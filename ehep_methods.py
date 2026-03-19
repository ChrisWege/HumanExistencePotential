"""
Reads in the position of archaeological sites, monthly temperature and
precipitation values and reads out the Environmental Human Existence Potential
computed by logistic regression.

Konstantin Klein
Institute of Geophysics and Meteorology
University of Cologne
konstantin.klein@uni-koeln.de

modifications:
2024-12-20  Annika Vogel    split into different files, here: only function with methods
"""
### Import Modules ###

from netCDF4 import Dataset
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import numpy.ma as ma
from geopy.distance import great_circle
import configure as cf
import ehep_util as eu
import ehep_inout as eio
import random
from datetime import datetime
from pathos.multiprocessing import ProcessingPool as Pool
from itertools import product
import sys


#########################################
##### FILL_LABELARRAYS ##################
# fill label-specific index,coordinate,predicor arrays
# - author: Christian Wegener, IGMK, UniKoeln
# - version: ???(from /data/hescor/owf/hep on 2024-12-04)
# - modification history:
#2024-01-17     Annika Vogel    move from pre_abs_sites as new function for all lables
###
def fill_labelarrays(lonlat_idx,lat,lon,data_set):

    # - initialize output arrays
    label_pred = []
    label_lat = np.zeros(len(lonlat_idx))
    label_lon = np.zeros(len(lonlat_idx))

    for i in range(len(lonlat_idx)):
        # - gridpoint indices
        y = lonlat_idx[i][0]
        x = lonlat_idx[i][1]

        # - cooridnates
        if cf.gridtype_t == 'curvilinear':
            label_lat[i] = lat[y,x]
            label_lon[i] = lon[y,x]
        else:
            label_lat[i] = lat[y]
            label_lon[i] = lon[x]

        # - predictor values
        if len(np.shape(data_set)) > 2:
            z = []
            for j in range(len(data_set)):
                z.append(data_set[j, y, x])
            label_pred.append(z)
        else:
            label_pred.append(data_set[y, x])

    return label_lat, label_lon, label_pred


#########################################
##### PRE_ABS_SITES #####################
# specify presence/pseudo-absence/apriori-absence points (training domain)
# - author: Christian Wegener, IGMK, UniKoeln
# - version: ???(from /data/hescor/owf/hep on 2024-12-04)
# - modification history:
#2024-12-20     Annika Vogel    change order of plotting presence/absense
#2025-01-17     Annika Vogel    generalize filling of label-specific arrays into function
#2025-02-27     Annika Vogel    check land-sea mask for domain
#2025-03-31     Annika Vogel    add plotting histogram
#2025-04-10     Annika Vogel    generalize labeling wrt missig bioclim fields
#2025-06-24     Annika Vogel    option for labeling apriori abs for loc exceeding limits of presence cond
###
def pre_abs_sites(lat, lon, dlat, dlon, data_set, sites, mainin_array, start):
    """
     Presence and pseudo-absence grid points will be calculated based on the
     archaeological sites. Grid points in specified radius around a site will
     be treated as presence, every other grid point, which is not farther away
     than 1000 km to the next site, is a pseudo-absence point
    :param sites position
    :return: 2 arrays with lat-lon values of presence (1) and
             pseudo-absence points (2)
     Stuff that needs to be changed in the absence calculations:
    #Load in Orography and check if elevation is below 1400m
    #If any of them gives a "false", the point is a-priori absence
    #In the calculations, we now have to distinguish these points:
    ### a-priori absence is always valid
    ### pseudo-absence is only valid in 1/3 of the times and is randomized
    ### within the 1000 repetitions, even before the 80/20 cut of the
    ### validation scheme.
    """
    print(' pre_radius_site=',cf.pre_radius_site)
    print(' absapri_limits_mode=',cf.absapri_limits_mode)
    print(' sample_factor=',cf.sample_factor)
    print(' train_absapri_only=',cf.train_absapri_only)
    print(' train_absapri_minnum=',cf.train_absapri_minnum)
    print(' dataratio_trai=',cf.dataratio_trai)
    print(' abs_fraction=',cf.abs_fraction)

    # - read land-sea mask from file
    data = Dataset(cf.path_land_sea_mask)
    print('Reading land-sea mask from file:',cf.path_land_sea_mask)
    data_lat = np.array(data.variables["lat"])
    data_lon = np.array(data.variables["lon"])
    land_sea_mask = np.array(data.variables["land_sea_mask"])
    data.close()
    #-check (training) domain
    eu.check_domain(cf.lat_min_t,cf.lat_max_t,cf.lon_min_t,cf.lon_max_t,data_lat,data_lon)

    # - initialize index arrays
    pres_indices = []
    apriabs_indices = []
    abs_indices = []

    #*** FUNCTION: label_gridpoints ***#
    # assign presence/pseudo-/apiori-absence/none label to single point
    def label_gridpoints(y,x):
        # - 3:none (sea,NaN)
        if eu.land_or_sea(data_lat,data_lon,land_sea_mask,lat[y],lon[x]) == 0:
            return [y,x,3]
        if np.isnan(mainin_array[0,y,x]): #non-normalized
            return [y,x,3]
        if np.isnan(mainin_array[-1,y,x]): #non-normalized
            return [y,x,3]

        # - 1:apiori-absence (exceed specified limits of bioclim fields)
        if cf.absapri_limits_mode == 1:
            for ibio in range(eu.nmainin):
                bionr_tmp = int(eu.mainin_fieldnr[ibio])
                if (bionr_tmp == 1 and \
                        ( mainin_array[ibio,y,x] > cf.bio1_max or mainin_array[ibio,y,x] < cf.bio1_min ) ): #non.normalized
                    #'bio1'
                    return [y,x,1]
                elif (bionr_tmp == 2 and \
                        ( mainin_array[ibio,y,x] > cf.bio2_max or mainin_array[ibio,y,x] < cf.bio2_min ) ): #non-normalized
                    #'bio2'
                    return [y,x,1]
                else:
                    continue

        # - 0:presence (within radius around site)
        for i in range(len(sites)):
            if (abs(lat[y] - sites[i, 0]) > (2*cf.sample_factor)) | (abs(lon[x] - sites[i, 1]) > (2*cf.sample_factor)):
                continue
            if great_circle((lat[y], lon[x]),
                                (sites[i, 0], sites[i, 1])).km <= cf.pre_radius_site:
                return [y,x,0]

        # - 2:pseudo-absence(else)
        return [y,x,2]
    #*** ***#

    #*** FUNCTION: relabel2apri_limits_any ***#
    # re-label points to apriori if exceeding given limits for ANY infield
    def relabel2apri_limits_any(y,x,label):
        if label != 3: #skip if sea/NaN
            #-check limits for every infield(normalized)
            for ifield in range(eu.ninfields_use):
                if data_set[ifield,y,x] < pres_min[ifield] or data_set[ifield,y,x] > pres_max[ifield]:
                    return [y,x,1] #change to apriori-absence
    
        #-if no infield exceeds limits
        return [y,x,label] #keep label
    #*** ***#

    #*** FUNCTION: relabel2apri_limits_all ***#
    # re-label points to apriori if exceeding given limits for ALL infields
    def relabel2apri_limits_all(y,x,label):
        if label == 3: #skip if sea/NaN
            return [y,x,label] #keep label

        #-check limits for every infield(normalized)
        for ifield in range(eu.ninfields_use):
            if data_set[ifield,y,x] >= pres_min[ifield] and data_set[ifield,y,x] <= pres_max[ifield]:
                return [y,x,label] #keep label(if any field lies within limits)

        #-if any infield exceeds limit
        return [y,x,1] #change to apriori-absence
    #*** ***#
    
    # - assign label to each point(parallel)
    latlon_yx = np.array([[y,x] for y in range(len(lat)) for x in range(len(lon))])
    pool = Pool(cf.process_count)
    gridlabels = pool.amap(label_gridpoints,latlon_yx[:,0],latlon_yx[:,1])
    gridlabels = np.array(gridlabels.get())

    #-re-label points to apriori if exceeding min/max of presence conditions
    if cf.absapri_limits_mode == 2 or cf.absapri_limits_mode == 3:
        print('re-labeling to apriori when exceeding limits of any presence cond...')
        #-(for re-label only!):extract gridpoint indices for presence
        for k in range(len(gridlabels)):
           if gridlabels[k,2] == 0:
                pres_indices.append([gridlabels[k,0],gridlabels[k,1]])

        #-(for re-label only!):extract coordinate and predictor array for presence
        pres_lat, pres_lon, pres_pred = fill_labelarrays(pres_indices,lat,lon,data_set)
        #-calc min/max at presence points for each infield
        pres_max = np.zeros(eu.ninfields_use)
        pres_min = np.zeros(eu.ninfields_use)
        for ifield in range(eu.ninfields_use):
            pres_ifield = [row[ifield] for row in pres_pred]
            pres_max[ifield] = np.nanmax(pres_ifield)
            pres_min[ifield] = np.nanmin(pres_ifield)
        #-re-label pseudo absence points(parallel)
        pool = Pool(cf.process_count)
        if cf.absapri_limits_mode == 2:
            gridlabels = pool.amap(relabel2apri_limits_any,gridlabels[:,0],gridlabels[:,1],gridlabels[:,2])
        else: #absapri_limits_mode == 3
            gridlabels = pool.amap(relabel2apri_limits_all,gridlabels[:,0],gridlabels[:,1],gridlabels[:,2])

        gridlabels = np.array(gridlabels.get())

    # - extract array with gridpoint indices for each label
    for k in range(len(gridlabels)):
        if gridlabels[k,2] == 0:
            pres_indices.append([gridlabels[k,0],gridlabels[k,1]])
        elif gridlabels[k,2] == 1:
            apriabs_indices.append([gridlabels[k,0],gridlabels[k,1]])
        elif gridlabels[k,2] == 2:
            abs_indices.append([gridlabels[k,0],gridlabels[k,1]])
        elif gridlabels[k,2] == 3:
            continue
        else:
            print("ERROR: Something went wrong in unpacking of the parallel script.")
            print(gridlabels[k])
            sys.exit()

    print("Presence, Absence and Apriori-Absence Points calculated")
    print(datetime.now() - start)
    
    # - extract coordinate and predictor array for each label
    pres_lat, pres_lon, pres_pred = fill_labelarrays(pres_indices,lat,lon,data_set)
    print("Presence/absence predictors calculated")
    print(datetime.now() - start)

    abs_lat, abs_lon, abs_pred = fill_labelarrays(abs_indices,lat,lon,data_set)
    apriabs_lat, apriabs_lon, apriabs_pred = fill_labelarrays(apriabs_indices,lat,lon,data_set)
    print("Predictors Absence calculated")
    print(datetime.now() - start)

    # --- plot presence-absence map
    if cf.plot_presabs:
        eio.plot_presabs(lat,lon,abs_lat,abs_lon,pres_lat,pres_lon,apriabs_lat,apriabs_lon,sites)

    # --- plot histogram of normalized fields
    if cf.plot_hist:
        eio.plot_histogram(pres_indices,abs_indices,apriabs_indices,lat,lon,mainin_array)


    return pres_pred, abs_pred, apriabs_pred, pres_indices, abs_indices, apriabs_indices


#########################################
##### TRAINING_TEST_SET #################
# prepare training and test data
# - author: Christian Wegener, IGMK, UniKoeln
# - version: ???(from /data/hescor/owf/hep on 2024-12-04)
# - modification history:
#2024-12-20     Annika Vogel    generalize: training/test ratio as config param
#                               limits for apriori absense points as config param
#2025-02-24     Annika Vogel    generalize for deterministic run (using all data for training)
#                               ->no test data, no statistical evaluation
#                               fraction of used speudo-absence data as config param
###

def training_test_set(data_set, pre_data, abs_data, apri_data, pre_points, abs_points, apri_points):
    """
    Use the data set and presence and absence locations to generate a training
    (defined amount of presence and pseudo-absence points) and test set (the
    rest of presence and pseudo-absence points)
    Use only one third of pseudo-absence points before taking the defined cut
    :return: training and test data
    """

    # - shuffle data randomly in array
    random.seed()
    predictors_presence = random.sample(pre_data,len(pre_data))
    pre_points_shuffled = random.sample(pre_points,len(pre_points))

    predictors_absence = random.sample(abs_data,len(abs_data))
    abs_points_shuffled = random.sample(abs_points,len(abs_points))
    ##-only use predefined fraction of speudo-absence points
    predictors_absence = predictors_absence[0:eu.nabs_use]
    abs_points_shuffled = abs_points_shuffled[0:eu.nabs_use]

    predictors_apriori_absence = random.sample(apri_data,len(apri_data))
    apri_points_shuffled = random.sample(apri_points,len(apri_points))

    # - check number of apriori absence points if train only on them
    if cf.train_absapri_only and (napriabs_trai < cf.train_absapri_minnum):
        print('CAUTION: number of apriori absence point too low -> use also pseudo absense for training!')
        train_absapri_only = False
    else:
        train_absapri_only = cf.train_absapri_only

    # - fill training data arrays
    #-get dimensions
    if cf.train_absapri_only:
        nabs_trai = 0
    else:
        nabs_trai = eu.nabs_trai

    ntrai = eu.npres_trai+nabs_trai+eu.napriabs_trai
    iapri_trai_start = eu.npres_trai+nabs_trai
    
    #-fill arrays
    if len(np.shape(data_set)) > 2:
        x_train = np.zeros([ntrai,len(data_set)])
        if eu.napriabs_trai > 0:
            x_train[iapri_trai_start:] = predictors_apriori_absence[0:eu.napriabs_trai]
        if not cf.train_absapri_only: #use also (pseudo) absence points
            x_train[eu.npres_trai:iapri_trai_start] = predictors_absence[0:eu.nabs_trai]

        x_train[0:eu.npres_trai] = predictors_presence[0:eu.npres_trai]

    else:
        x_train = np.zeros([ntrai,1])
        if eu.napriabs_trai > 0:
            x_train[iapri_trai_start:,0] = predictors_apriori_absence[0:eu.napriabs_trai]
        if not cf.train_absapri_only: #use also (pseudo) absence points
            x_train[eu.npres_trai:iapri_trai_start, 0] = predictors_absence[0:eu.nabs_trai]

        x_train[0:eu.npres_trai, 0] = predictors_presence[0:eu.npres_trai]

    y_train = np.zeros(ntrai)
    y_train[0:eu.npres_trai] += 1

    # - fill testing data arrays
    #-get dimensions
    if cf.train_absapri_only:
        nabs_test = 0
    else:
        nabs_test = eu.nabs_test

    ntest = eu.npres_test+nabs_test+eu.napriabs_test
    iapri_test_start = eu.npres_test+nabs_test
    #-fill arrays
    if len(np.shape(data_set)) > 2:
        x_test = np.zeros([ntest,len(data_set)])
        if eu.napriabs_test > 0:
            x_test[iapri_test_start:] = predictors_apriori_absence[eu.napriabs_trai:]
        if not cf.train_absapri_only: #use also (pseudo) absence points
            if eu.nabs_test > 0:
                x_test[eu.npres_test:iapri_test_start] = predictors_absence[eu.nabs_trai:]

        if eu.npres_test > 0:
            x_test[0:eu.npres_test] = predictors_presence[eu.npres_trai:]

    else:
        x_test = np.zeros([ntest,1])
        if eu.napriabs_test > 0:
            x_test[iapri_test_start:,0] = predictors_apriori_absence[eu.napriabs_trai:]
        if not cf.train_absapri_only: #use also (pseudo) absence points
            if eu.nabs_test > 0:
                x_test[eu.npres_test:iapri_test_start, 0] = predictors_absence[eu.nabs_trai:]

        if eu.npres_test > 0:
            x_test[0:eu.npres_test, 0] = predictors_presence[eu.npres_trai:]

    y_test = np.zeros(ntest)
    if eu.npres_test > 0:
        y_test[0:eu.npres_test] += 1

    # - define matrix of presence/absence indices in domain
    if len(np.shape(data_set)) > 2:
        pre_abs_matrix = ma.zeros([np.shape(data_set)[1], np.shape(data_set)[2]])
    else:
        pre_abs_matrix = ma.zeros([np.shape(data_set)])

    pre_abs_matrix[:] = ma.masked

    for k in range(eu.npres_trai):
        pre_abs_matrix[pre_points_shuffled[k][0], pre_points_shuffled[k][1]] = 1
    for k in range(eu.nabs_trai):
        pre_abs_matrix[abs_points_shuffled[k][0], abs_points_shuffled[k][1]] = 0
    if eu.napriabs_trai > 0:
        for k in range(eu.napriabs_trai):
            pre_abs_matrix[apri_points_shuffled[k][0], apri_points_shuffled[k][1]] = 0

    #print("NAN? ",np.count_nonzero(np.isnan(x_train)))
    x_train[np.isnan(x_train)] = 0.
    #print("NAN? ",np.count_nonzero(np.isnan(x_train)))
    y_train[np.isnan(y_train)] = 0.
    x_test[np.isnan(x_test)] = 0.
    y_test[np.isnan(y_test)] = 0.

    return x_train, y_train, x_test, y_test, pre_abs_matrix


#########################################
##### TRAIN_LOGREG ######################
# train logistic regression to data
# - author: Christian Wegener, IGMK, UniKoeln
# - version: ???(from /data/hescor/owf/hep on 2024-12-04)
# - modification history:
# 2025-03-04    Annika Vogel    initialize polynomials only once, get poly names
###
def train_logreg(x_train, y_train, x_test, y_true, iv_set, dlat, dlon, hep_mask):
    """
    Training the logistic regression model with train dataset features (train_x)
    and target (train_y)
    """
    # - initialize 2nd degree polynomials
    poly = PolynomialFeatures(degree=2)

    # - perform logistic regression for 2nd degree polynomials
    #-2nd order polynomials of each training input vector (=diff predictors at each location)
    x_train = eu.second_degree_polynomials(poly,x_train)
    eu.trainpoly_names = poly.get_feature_names_out(input_features=eu.trainfield_names)
    trained_logreg = LogisticRegression(penalty='l1', C=cf.logreg_lasso, fit_intercept=False,
                                        solver='saga', tol=cf.logreg_tol, max_iter=cf.logreg_max_iter,
                                        class_weight='balanced')
    trained_logreg.fit(x_train, y_train)
    # - get coefficients of each polynomial term
    coefs = trained_logreg.fit(x_train, y_train).coef_[0]
    coefs = np.array(coefs)

    # - calc modelled outcome (here: HEP) for (2nd order poly of) investigation data
    iv_set = eu.second_degree_polynomials(poly,iv_set)
    #-get probability outcome
    prediction = trained_logreg.predict_proba(iv_set)
    phi_e = eu.back_transpose_data_set(prediction[:, 1], dlat, dlon)
    phi_e.mask = hep_mask
    phi_e.filled(0.)

    # - statistical evaluation of model with (2nd order poly of) testing data
    if x_test.shape[0] > 0:
        x_test = eu.second_degree_polynomials(poly,x_test)
        #-get probability outcome
        y_test = trained_logreg.predict_proba(x_test)[:, 1]
        #-brier score
        y_test_brier = 1. / (1 + np.exp(-coefs[0] * x_test[:, 0]))
        brier_score1 = sum((y_true - y_test) ** 2)
        brier_score2 = sum((y_true - y_test_brier) ** 2)
        brier_score = 1 - brier_score1 / brier_score2
        #-area-under-curve
        auc = roc_auc_score(y_true, y_test)
    else:
        brier_score = -1.
        auc = -1.

    return phi_e, coefs, auc, brier_score


#########################################
##### TRAIN_RANDOMFOREST ################
# train random forest to data
# - author: Christian Wegener, IGMK, UniKoeln
# - version: ???(from /data/hescor/owf/hep on 2024-12-04)
# - modification history:
#2025-01-21     Annika Vogel    rebugging handling of RandomForest
###
def train_RandomForest(x_train, y_train, x_test, y_true, iv_set, dlat, dlon, hep_mask):
    """
    Training the random forest model with train dataset features (train_x)
    and target (train_y)
    :return: phi_e, feature_importance, auc, brier_score, oob_score
    """

    # - train random forest model
    rf = RandomForestClassifier(n_estimators=150, random_state=0,
                                oob_score=True, class_weight='balanced')
    rf.fit(x_train, y_train)

    # - calc modelled outcome (here: HEP) for investigation data
    #-get probability outcome
    prediction = rf.predict_proba(iv_set)
    phi_e = eu.back_transpose_data_set(prediction[:, 1], dlat, dlon)
    phi_e.mask = hep_mask
    phi_e.filled(0.)
    features = rf.feature_importances_

    # - statistical evaluation of model with testing data
    #-get probability outcome
    y_test = rf.predict_proba(x_test)[:, 1]
    #-brier score
    brier_score = sum((y_true - y_test) ** 2) / len(y_true)
    #-area-under-curve
    auc = roc_auc_score(y_true, y_test)
    oob_score = rf.oob_score_

    return phi_e, features, auc, brier_score, oob_score



#########################################
##### FIT_SIMPLE_GAUSSIAN ################
# fit simple Function (Gaussian) to data
#  -mean, std of each field over presence points
#  -total HEP=weighted multipl.mean(Gaussians): HEP(i,j) = exp[1/N*sum_n( -(x_n(i,j)-mean_n)^2 / 2*std_n^2 )] ,n=[1,N]:input fields
# - author: Annika Vogel Wegener, IGMK, UniKoeln
# - modification history:
#2025-03-06     Annika Vogel    initial version
#2025-04-22     Annika Vogel    calc distinctiveness, use as weight for input fields
#2025-04-23     Annika Vogel    add calc of skewness (for Gamma distr)
###
def train_simple_fit(pres,iv_set,hep_mask,dlat,dlon):
    # - load modules
    from scipy.stats import skew

    # - ini prediction fields
    npoints = iv_set.shape[0]
    pres_means = np.zeros(eu.ninfields_use)
    pres_stds = np.zeros(eu.ninfields_use)
    pres_skews = np.zeros(eu.ninfields_use)
    pres_gamma_alpha = np.zeros(eu.ninfields_use)
    pres_gamma_beta = np.zeros(eu.ninfields_use)
    pres_gamma_shift = np.zeros(eu.ninfields_use)
    lnphi_i = np.zeros(iv_set.shape) #(npoints,ninfields_use)
    weight = np.zeros(eu.ninfields_use)
    phi_e = np.zeros(npoints)

    # - ini distinctiveness fields
    dmeans = np.zeros(eu.ninfields_use)
    dstds = np.zeros(eu.ninfields_use)
    distinct = np.zeros(eu.ninfields_use)

    # --- calc statistics at presence points, distinctiveness (all .vs. pres) of human presence conditions for each input field

    for iuse in range(eu.ninfields_use):

        # - get array
        pres_column = [row[iuse] for row in pres]

        # - calc stats
        pres_means[iuse] = np.nanmean(pres_column)
        pres_stds[iuse] = np.nanstd(pres_column)
        pres_skews[iuse] = skew(pres_column, nan_policy='omit')
        # - calc param of Gamma distr
        pres_gamma_alpha[iuse] = 4./pres_skews[iuse]**2
        pres_gamma_beta[iuse] = 0.5*pres_stds[iuse]*pres_skews[iuse]
        pres_gamma_shift[iuse] = pres_means[iuse] - 2.*pres_stds[iuse]/pres_skews[iuse]
        print('normalized ',eu.trainfield_names[iuse],' @presence points: mean=',pres_means[iuse],', std=',pres_stds[iuse],', skew=',pres_skews[iuse], \
                ',-> a=',pres_gamma_alpha[iuse],', b=',pres_gamma_beta[iuse],', c=',pres_gamma_shift[iuse])

        # - calc distinctiveness
        #-difference in means (mean over all points in training domain =0 per def of normalization)
        dmeans[iuse] = pres_means[iuse] - 0.
        #-difference in standard deviations (std over all points in training domain =1 per def of normalization)
        dstds[iuse] = pres_stds[iuse] - 1.
        #-total distinctiveness: root mean square(delta mean,delta stdev)
        #distinct[iuse] = dmeans[iuse]**2+dstds[iuse]**2 #np.sqrt
        #(alternative: use only stdev->variance, ignore all fields with stdev>0)
        distinct[iuse] = dmeans[iuse]**2+min(dstds[iuse],0.)**2 #=...+max(-dstds,0)**2
        #distinct[iuse] = max(-dstds[iuse],0.)**2 #=only use negative dstdev =smaller stdev for pres
        if dstds[iuse] >= 0.:
            print('normalized stdev of ',eu.trainfield_names[iuse],' larger for presence-points than for total domain ->reduced weight!')

    distinct_sum = np.sum(distinct)
    for iuse in range(eu.ninfields_use):
        # - calc normalized weight wrt input fields
        #weight[iuse] = 1./eu.ninfields_use #(org,no weight)
        weight[iuse] = distinct[iuse]/distinct_sum
        # - application to investigation data: HEP component for each field
        lnphi_i[:,iuse] = - weight[iuse]* ((iv_set[:,iuse]-pres_means[iuse])**2) / (2*(pres_stds[iuse]**2))

        print('distincit(mean-part)=',dmeans[iuse]**2,', std-part=',min(dstds[iuse],0.)**2)
        print(eu.trainfield_names[iuse],': distinctiveness(rms)=',distinct[iuse], ' ->weight=',weight[iuse])
        #print('input field weights=',weight)

    # - finalize total prediction
    phi_1d = np.exp(np.sum(lnphi_i,axis=1)) #(npoints)
    phi_e = eu.back_transpose_data_set(phi_1d[:], dlat, dlon) #(nx,ny)
    #prediction[:, 1]
    phi_e.mask = hep_mask #(nx,ny)
    phi_e.filled(0.)
    print('HEP range: min=',np.min(phi_e),', mean=',np.mean(phi_e),', max=',np.max(phi_e))

    # --- plot distinctiveness
    if cf.plot_distinct: eio.plot_distinct(pres_means,pres_stds,distinct)

    return phi_e


#########################################
##### MAIN_CALCULATION ##################
# do main calculation of HEP (parallel)
# - author: Christian Wegener, IGMK, UniKoeln
# - version: ???(from /data/hescor/owf/hep on 2024-12-04)
# - modification history:
#2025-01-21     Annika Vogel    rebugging handling of RandomForest
###
def main_calculation(irun, iv_set, dlat, dlon, training_set, pre, abs, apri, pre_points, abs_points, apri_points,hep_mask):
    """

    :return:
    """

    # - select training/test data (random)
    x_train, y_train, x_test, y_true, pre_abs_matrix = \
        training_test_set(training_set, pre, abs, apri, pre_points, abs_points, apri_points)

    # - train model fit on training data
    if cf.model_training == 'rf':
        #-random forest
        phi_e, features, auc, brier_score, obb_score = \
            train_RandomForest(x_train, y_train, x_test, y_true, iv_set, dlat, dlon, hep_mask)
        return phi_e, features, auc, brier_score, pre_abs_matrix #, obb_score

    elif cf.model_training == 'simpleFit':
        #-simple Gaussian fit
        phi_e = train_simple_fit(pre, iv_set, hep_mask, dlat, dlon)

        return phi_e, [0.], 0., 0., pre_abs_matrix

    else: #'logreg'
        #-logistic regression
        phi_e, coef, auc, brier_score = \
            train_logreg(x_train, y_train, x_test, y_true, iv_set, dlat, dlon, hep_mask)
        if irun == 0:
            print('2nd order polynomials: ',eu.trainpoly_names)

        return phi_e, coef, auc, brier_score, pre_abs_matrix

