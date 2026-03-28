import numpy as np
import numpy.ma as ma
import configure as cf

"""
Provides initialization of global dimensions and statistics, data
transpose/back-transpose helpers, domain and land/sea checks, and
polynomial feature construction used across the HEP model.
"""

##### INI_GLOBAL_PRE ####################
# initialize global variables for reading (to be used in all functions across files)
# !BEFORE reading: temporal dim/arrays, to be reduced later to available fields-only!

def ini_global_pre():
    global runs, nmainin, ninfields, ninfields_use
    global idx_use_in_all
    global mainin_fieldnr, mainin_fieldnr_use

    # - check for valid input parameters: number of parallel runs
    if cf.dataratio_trai == 1. and cf.abs_fraction == 1.:
        #-no stochastic selection, reduce to single deterministic run
        print('NOTE: no stochastic selection ->Only perform sinlge deterministic run!')
        runs = 1
    else:
        # - set of stochastic runs
        runs = cf.runs

    # - PRE-set dimension of bioclim fields(to be updated after reading fields [ini_global_post]
    # if cf.input_filedim_type_t == 'time':
    #    nmainin = len(cf.input_filetime_fieldnr)
    if cf.input_onefield_t:
        nmainin = cf.nbioclim_def
    else:
        nmainin = len(cf.input_varnames_t)

    # - PRE-set dimension of all input fields (to be updated after reading fields)
    ninfields = nmainin
    ninfields_use = len(cf.input_var_use)
    if cf.soil_use:
        ninfields = ninfields+1
        ninfields_use = ninfields_use+1

    # - PRE-set list of available field numbers of main input (starting with 1!)
    mainin_fieldnr = np.zeros(nmainin)

    # - PRE-set array with used main input fieldnames
    mainin_fieldnr_use = np.zeros(len(cf.input_var_use))

    # - PRE-set list of indices of used input fields in all input fields
    idx_use_in_all = np.zeros(ninfields_use)


##### INI_GLOBAL_POST ###################
# initialize global variables for calculation (to be used in all functions across files)
# !AFTER reading: reduce existing dim/arrays wrt available fields, and ini all remaining!

def ini_global_post(mainidx_to_rm,mainin_array):
    global nmainin, ninfields
    global mainin_fieldnr
    global allfield_names
    global training_mean_all, training_stdev_all

    # - reduce existing dimensions and arrays
    if mainidx_to_rm:
        nmainin -= len(mainidx_to_rm)
        ninfields -= len(mainidx_to_rm)
        for imain in mainidx_to_rm:
            np.delete(mainin_array,imain,0) #delete imain-th row from 3D array

    # - initialize names of training fields
    allfield_names = ['']*ninfields

    # - initialze normalization fields
    training_mean_all = np.zeros(ninfields)
    training_stdev_all = np.zeros(ninfields)

    return mainin_array


##### INI_GLOBAL_USE ####################
# initialize global variables related to _use for calculation (to be used in all functions across files)

def ini_global_use():
    global ninfields_use, ninfields_poly
    global trainfield_names, trainpoly_names
    global ninfields_useext
    global training_mean_use, training_stdev_use

    # - initialize other dimensions
    ninfields_poly = int((ninfields_use*(ninfields_use+1))/2 +ninfields_use +1)
    print(' eu.ninfields_poly=',ninfields_poly)
    # - if simpleFit and infields_ext_mode==1: ini extended used fields
    if cf.model_training == 'simpleFit' and cf.infields_ext_mode == 1:
        ninfields_useext = int( (ninfields_use*(ninfields_use+1))/2 )
    else:
        ninfields_useext = ninfields_use

    # - initialize names and stats of training fields
    trainpoly_names = ['']*ninfields_poly
    trainfield_names = ['']*ninfields_useext
    training_mean_use = np.zeros(ninfields_useext)
    training_stdev_use = np.zeros(ninfields_useext)


##### INI_DIMS_USE ######################
# initialize global dimensions of used training data

def ini_dims_use(npres, nabs, napriabs):
    global abs_fraction, dataratio_trai
    global npres_use, nabs_use, napriabs_use
    global npres_trai, nabs_trai, napriabs_trai
    global npres_test, nabs_test, napriabs_test

    # - check for valid input parameters
    # - fraction of used pseudo-absence points
    if cf.abs_fraction <= 0. or cf.abs_fraction > 1.:
        abs_fraction = 1./3.
        print('NOTE: invalid fraction of used pseudo-absence points. Set to 1/3')
    else:
        abs_fraction = cf.abs_fraction
    # - training/test ratio
    if cf.dataratio_trai <= 0. or cf.dataratio_trai > 1.:
        dataratio_trai = 0.8
        print('NOTE: invalid dataratio_trai and/or dataratio_test. Set to 0.8')
    else:
        dataratio_trai = cf.dataratio_trai

    # - get number of used data
    npres_use = npres
    napriabs_use = napriabs
    # - only use predefined fraction of speudo-absence points
    nabs_use = int(np.floor(nabs*abs_fraction))

    # - get number of training/test data
    npres_trai = int(np.ceil(npres_use * dataratio_trai))
    nabs_trai = int(np.ceil(nabs_use * dataratio_trai))
    napriabs_trai = int(np.ceil(napriabs_use * dataratio_trai))
    npres_test = npres_use - npres_trai
    nabs_test = nabs_use - nabs_trai
    napriabs_test = napriabs_use - napriabs_trai



##### TRANSPOSE_DATA_SET ################
# transpose data from 3D(var,lat,lon) to 2D(lat*lon,var)
# !CAUTION: lat, lon must be in 2nd and 3rd dimension of data_SET!

def transpose_data_set(data_set, dlat, dlon):
    """
    Transpose the data set to fit the predict tool of the logistic regression:
    """
    if len(np.shape(data_set)) > 2:
        nvars = data_set.shape[0] #length of 1st dimension
        data_set_transposed = ma.zeros([dlat * dlon, nvars])
        for i in range(dlon):
            for j in range(nvars):
                data_set_transposed[range(i * dlat, (i + 1) * dlat), j] = data_set[j, :, i]
    else:
        data_set_transposed = ma.zeros([dlat * dlon, 1])
        for i in range(dlon):
            data_set_transposed[range(i * dlat, (i + 1) * dlat), 0] = data_set[:, i]
    return data_set_transposed


##### BACK_TRANSPOSE_DATA_SET ###########
# back-transpose data from (lat*lon,:) to 3D(:,lat,lon)

def back_transpose_data_set(data_set, dlat, dlon):
    """
    Back transpose the predicted data set in lat-lon form
    """
    back_transposed = ma.zeros([dlat, dlon])
    for i in range(dlon):
        back_transposed[:, i] = data_set[range(i * dlat, (i + 1) * dlat)]

    return back_transposed


##### CHECK_DOMAIN ########################
# check if lat/lon fields cover defined domain

def check_domain(lat_min,lat_max,lon_min,lon_max,lat_set,lon_set):
    """
    :return:
    """
    ### CHECK DOMAIN ###
    # - get min and max lat/lon: rounded to +-0.1deg
    lat_set_min = ma.floor(np.min(lat_set)*10)/10
    lat_set_max = ma.ceil(np.max(lat_set)*10)/10
    lon_set_min = ma.floor(np.min(lon_set)*10)/10
    lon_set_max = ma.ceil(np.max(lon_set)*10)/10

    # - check if completely outside domain
    if (lat_set_min >= lat_max) | (lat_set_max <= lat_min) | \
            (lon_set_min >= lon_max) | (lon_set_max <= lon_min):
        print("ERROR: Chosen domain non compatible with the"
              " data set.")
        print('lat_set_min(=',lat_set_min,') ?>? lat_max(=',lat_max,')')
        print('lat_set_max(=',lat_set_max,') ?<? lat_min(=',lat_min,')')
        print('lon_set_min(=',lon_set_min,') ?>? lon_max(=',lon_max,')')
        print('lon_set_max(=',lon_set_max,') ?<? lon_min(=',lon_min,')')
        sys.exit()
    # - check if not fully inside domain
    elif (lat_set_min > lat_min) | (lat_set_max < lat_max) | \
            (lon_set_min > lon_min) | (lon_set_max < lon_max):
        print("Warning: Data set not complete in chosen domain.")
        print('lat_set_min(=',lat_set_min,') ?<? lat_min(=',lat_min,')')
        print('lat_set_max(=',lat_set_max,') ?>? lat_max(=',lat_max,')')
        print('lon_set_min(=',lon_set_min,') ?<? lon_min(=',lon_min,')')
        print('lon_set_max(=',lon_set_max,') ?>? lon_max(=',lon_max,')')


##### LAND_OR_SEA #######################
# index of land (=1) or sea (=0) for individual lat-lonn coordinate

def land_or_sea(data_lat, data_lon, lsm, llat, llon):

    min_lat = np.where(abs(data_lat - llat) == np.min(abs(data_lat - llat)))[0][0]
    min_lon = np.where(abs(data_lon - llon) == np.min(abs(data_lon - llon)))[0][0]

    if lsm[min_lat, min_lon] == 1:
        return 1
    else:
        return 0


##### CUT_UNSPECIFIED_REGIONS ###########
# (currently unused!)
def cut_unspecified_regions(llat, llon):

    unspecified_point = 0
    if llat > 65.:
        unspecified_point += 1
        return unspecified_point

    if llon < -20.:
        unspecified_point += 1
    elif llon < 0:
        if llat < 35.95:
            unspecified_point += 1
    elif llon < 1.7:
        if llat < 37:
            unspecified_point += 1
    elif llon < 8:
        if llat < 40:
            unspecified_point += 1
    elif llon < 45.:
        if llat < 46:
            unspecified_point += 1
    else:
        unspecified_point += 1
    return unspecified_point


##### AFRICA_OR_EUROPE ##################
# (currently unused!)
def Africa_or_Europe(llat, llon):
    """
    Is the lattice point in Africa or Europe?
    :return:Africa = 1, Europe = 0
    """
    continent = 0

    if llon < 0:
        if llat < 35.95:
            continent += 1
    elif llon < 11.5:
        if llat < 37.5:
            continent += 1
    elif llon < 30:
        if llat < 33.5:
            continent += 1

    return continent


##### SECOND_DEGREE_POLYNOMIALS #########
# create 2nd degree polynomial of input vector
#     Example: input: array = ([a,b,c],[d,e,f])
#             output: poly_array = ([1,a,b,c,a^2,a*b,a*c,b^2,b*c,c^2],
#                                   [1,d,e,f,d^2,d*e,d*f,e^2,e*f,f^2])

def second_degree_polynomials(poly,data_set):

    polynomials = poly.fit_transform(data_set)

    return polynomials

