"""
# - author: Christian Wegener, IGMK, UniKoeln
# - version: ???(from /data/hescor/owf/hep on 2024-12-04)
# - modification history:
#2024-12-20     Annika Vogel    move functions from ehep_run.py in new file
"""

import configure as cf
import ehep_util as eu
import sys
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
#########################################
##### READ_DATA_NCDF ####################
# getting training/investigation data from netCDF file
# - author: Christian Wegener, IGMK, UniKoeln
# - version: ???(from /data/hescor/owf/hep on 2024-12-04)
# - modification history:
#2024-12-20     Annika Vogel    export into new function, call for training and investigation
#2025-01-13     Annika Vogel    only do spatial subsampling if needed
#2025-01-15     Annika Vogel    get number of input fields as input parameter
#2025-02-24     Annika Vogel    generalize: optional all bioclim variables in one field
#2025-02-25     Annika Vogel    generalize domain selection for different lat/lon directions
#2025-02-28     Annika Vogel    move domain check to new function, global def of normalization from training
#2025-03-04     Annika Vogel    define and print out bioclim names, print normalization from training
#2025-03-31     Annika Vogel    generalize for one bioclim field from different times in input file 'input_filedim_type=time'
#2025-04-09     Annika Vogel    generalize for any missing bioclim infput file
#2025-06-26     Annika Vogel    include option to extend input fields to cross-terms (simpleFit-only)
###
def read_data_ncdf(use_type,gridtype,lat_min,lat_max,lon_min,lon_max,input_path,input_latname,input_lonname,input_varnames,input_filedim_type,input_onefield,soil_path,soil_varname):
    """

    :return:
    """

    ### GET INPUT FIELDS ###
    # - get lat-lon fields and dimensions from main input file
    if input_filedim_type == 'time':
        input_path_tmp = input_path+cf.input_filetime_pathfield[0]+cf.input_filetime_pathend
        input_data = Dataset(input_path_tmp)
        #-prepare time index in file for reading main input fields at selected time
        times_set = input_data.variables[cf.input_filetime_timename]
        timeindex = -1
        if isinstance(times_set,str):
            #-if time indentifier is string
            for timeindex_tmp, time_tmp in enumerate(times_set):
                if cf.input_filetime_timeid in time_tmp:
                    timeindex = timeindex_tmp
                    break
        else:
            for timeindex_tmp, time_tmp in enumerate(times_set):
                if cf.input_filetime_timeid == time_tmp:
                    timeindex = timeindex_tmp
                    break
        print('used time for input file: timeindex=',timeindex,', times_set=',times_set[timeindex])
    else: #'field'
        input_data = Dataset(input_path)

    lat_set = input_data.variables[input_latname][:]
    lon_set = input_data.variables[input_lonname][:]
    #-adjust if longitude in [0;360]
    if np.max(lon_set) > 180.:
        lon_set = lon_set -180.

    print('domain lat:',lat_min,'°N - ',lat_max,'°N , lon:',lon_min,'°E - ',lon_max,'°E')

    '''
    # - CAUTION TEST-only:test 'curvilinear' with 'latlon' input...
    print("!!!CAUTION: USING 'latlon' DATA TO CHECK 'curvilinear' CODE!!!")
    dlat_set = len(lat_set)
    dlon_set = len(lon_set)
    lat_set_tmp = np.zeros([dlat_set,dlon_set])
    lon_set_tmp = np.zeros([dlat_set,dlon_set])
    for y in range(dlat_set):
        lat_set_tmp[y,:] = lat_set[y]
    for x in range(dlon_set):
        lon_set_tmp[:,x] = lon_set[x]
    lat_set = lat_set_tmp
    lon_set = lon_set_tmp
    #...TEST
    '''

    # - check if defined domain is within field domain
    if gridtype == 'curvilinear':
        dlat_set = np.shape(lat_set)[0]
        dlon_set = np.shape(lat_set)[1]
    else:
        dlat_set = len(lat_set)
        dlon_set = len(lon_set)

    print("read latitude points: ",dlat_set)
    print("read longitude points: ",dlon_set)
    eu.check_domain(lat_min,lat_max,lon_min,lon_max,lat_set,lon_set)

    input_data.close()

    # - select coordiantes within domain
    # - get indices of domain boundaries in fields
    # - update dimensions
    if gridtype == 'curvilinear':   #(CAUTION: to be tested!)
        # - (lat,lon) fields of each point in domain (wrt each point ->not yet implemented below!)
        #lat_domain = lat_set[(lat_set >= lat_min) & (lat_set <= lat_max) & (lon_set >= lon_min) & (lon_set <= lon_max)]
        #lon_domain = lon_set[(lat_set >= lat_min) & (lat_set <= lat_max) & (lon_set >= lon_min) & (lon_set <= lon_max)]

        #####
        #-extract indices of potential rows/columns in domain
        ny = 0 #ini
        y_domain = np.zeros(dlat_set)
        for y in range(dlat_set):
            if (np.max(lat_set[y,:]) >= lat_min) & (np.min(lat_set[y,:]) <= lat_max):
                y_domain[ny] = y
                ny += 1
        
        nx = 0 #ini
        x_domain = np.zeros(dlon_set)
        for x in range(dlon_set):
            if (np.max(lon_set[:,x]) >= lon_min) & (np.min(lon_set[:,x]) <= lon_max):
                x_domain[nx] = x
                nx += 1

        y_domain = y_domain[0:ny]
        x_domain = x_domain[0:nx]

        #-get domain boundary indices (CAUTION: to be tested!)
        #(assumes sorted coordinated, any order)
        y_0 = int(np.min(y_domain))
        y_end = int(np.max(y_domain))
        x_0 = int(np.min(x_domain))
        x_end = int(np.max(x_domain))

        #-reduce coordinate fields and update dimensions (CAUTION: to be tested!)
        lat_set = lat_set[y_0:y_end+1,x_0:x_end+1]
        lon_set = lon_set[y_0:y_end+1,x_0:x_end+1]
        #-2D dimensions of pther fields (:,dlat_set,dlon_set)
        dlat_set = ny
        dlon_set = nx
        #-dimension of lat/lon field (dlat_set_out)/(dlon_set_out)
        dlat_set_out = lat_set.shape[0]
        dlon_set_out = lat_set.shape[1]
        ######

        '''
        #-original(outdated!)
        y_0 = 0
        y_end = dlat_set - 1
        x_0 = 0
        x_end = dlon_set - 1
        c = 0
        for y in range(dlat_set):
            if np.max(lat_set[y, :]) >= lat_min:
                if c == 0:
                    c += 1
                    y_0 = y
            if np.min(lat_set[y, :]) >= lat_max:
                y_end = y
                break
        c = 0
        for x in range(dlon_set):
            if np.max(lon_set[:, x]) >= lon_min:
                if c == 0:
                    c += 1
                    x_0 = x
            if np.min(lon_set[:, x]) >= lon_max:
                x_end = x
                break

        lat_set = lat_set[y_0:y_end+1,x_0:x_end+1]
        lon_set = lon_set[y_0:y_end+1,x_0:x_end+1]
        dlat_set = np.shape(lat_set)[0]
        dlon_set = np.shape(lat_set)[1]
        '''

    else:
        #-(lat), (lon) fields within domain
        lat_y = lat_set[(lat_set >= lat_min) & (lat_set <= lat_max)]
        lon_x = lon_set[(lon_set >= lon_min) & (lon_set <= lon_max)]

        #-domain boundaries
        y_0 = np.where(lat_set==lat_y[0])[0][0]
        y_end = np.where(lat_set==lat_y[-1])[0][0]
        x_0 = np.where(lon_set==lon_x[0])[0][0]
        x_end = np.where(lon_set==lon_x[-1])[0][0]
        lat_set = lat_y #=lat_set[y_0:y_end+1]
        lon_set = lon_x #=lon_set[x_0:x_end+1]
        #-2D dimensions of other fields (:,dlat_set,dlon_set)
        dlat_set = len(lat_set)
        dlon_set = len(lon_set)
        #-dimension of each lat/lon field (dlat_set_out)/(dlon_set_out)
        dlat_set_out = dlat_set
        dlon_set_out = dlon_set

    print("use latitude points: ",dlat_set)
    print("use longitude points: ",dlon_set)


    # - read all input fields for selected domain from file
    # - attach other fields to single input array
    mainin_array = ma.zeros([eu.ninfields, dlat_set, dlon_set])

    # - main input (bioclim/vegetation)
    imain_avail = 0
    mainidx_to_rm = [] #list of main infield indices to be removed from dimensions and arrays (not available)
    if input_filedim_type == 'time':
        # - different times in one file, only one field in file (lat,lon)
        for imain in range(eu.nmainin):
            #-create field-specific filename
            input_path_itime = input_path+cf.input_filetime_pathfield[imain]+cf.input_filetime_pathend
            try:
                input_data_itime = Dataset(input_path_itime)
                mainin_array[imain_avail] = input_data_itime.variables[input_varnames[imain]][timeindex,y_0:y_end+1, x_0:x_end+1]
                input_data_itime.close()
            except:
                print('!!!CAUTION: main input file is not avialable ',input_path_itime)
                mainidx_to_rm.append(imain)
            else:
                eu.mainin_fieldnr[imain_avail] = int(imain+1)
                print('  -reading available input field ',input_varnames[imain],' from file ',input_path_itime)
                imain_avail += 1
                

    else: #'field'
        # - all fields in one file, only one time in file
        try:
            input_data = Dataset(input_path)
        except:
            print('!!!CAUTION: main input file is not avialable ',input_path)

        if input_onefield or len(input_varnames) == 1:
            #-all input variables in one field(var,lat,lon)
            for imain in range(eu.nmainin):
                try:
                    mainin_array[imain_avail] = input_data.variables[input_varnames[0]][imain,y_0:y_end+1, x_0:x_end+1]
                except:
                    print('!!!CAUTION: main input field is not aviailable, field=',input_varnames[0],'[',imain,',:,:]')
                    mainidx_to_rm.append(imain)
                else:
                    eu.mainin_fieldnr[imain_avail] = int(imain+1)
                    print('  -reading available input field ',input_varnames[0],'[',imain,',:,:]')
                    imain_avail += 1

        else:
            #-each variable in separate field(lat,lon)
            for imain in range(eu.nmainin):
                try:
                    mainin_array[imain_avail] = input_data.variables[input_varnames[imain]][y_0:y_end+1, x_0:x_end+1]
                except:
                    print('!!!CAUTION: input field is not aviailable, field=',input_varnames[imain],'[:,:]')
                    mainidx_to_rm.append(imain)
                else:
                    eu.mainin_fieldnr[imain_avail] = int(imain+1)
                    print('  -reading avialable input field ',input_varnames[imain])
                    imain_avail += 1

        input_data.close()

    # - reduce dimensions and arrays to available input fields(for all filedim_type)
    if use_type == 't': mainin_array = eu.ini_global_post(mainidx_to_rm,mainin_array)

    # - other (here: soil)
    if cf.soil_use:
        print('Reading soil data from file(s):',soil_path)
        isoilarray = eu.nmainin #updated(not available input fields removed)
        #eu.idx_use_in_all[imain_use] = isoilarray
        soil_data = Dataset(soil_path)
        soil_array = soil_data.variables[soil_varname][y_0:y_end+1, x_0:x_end+1]
        soil_data.close()
        mainin_array[isoilarray] = soil_array

    #mainin_array[np.where(mainin_array > 1e+30)] = ma.masked
    mainin_array = ma.masked_invalid(mainin_array)
    mainin_array.filled(0.)
    #if _t:
    mainin_array_full = ma.copy(mainin_array)
    #if _i:
    mask = mainin_array[0].mask

    # - only for training data: define field names and normalization values
    if use_type == 't':
        #-get field names
        for i in range(eu.ninfields):
            if cf.soil_use and i == isoilarray:
                eu.allfield_names[i] = cf.soil_varname_t
            elif cf.input_onefield_t:
                eu.allfield_names[i] = 'bio'+str(int(eu.mainin_fieldnr[i])) #str(i+1)
            else:
                eu.allfield_names[i] = cf.input_varnames_t[i]

        #-calc normalization: domain mean and std of each field
        for i in range(eu.ninfields):
            eu.training_mean_all[i] = np.mean(mainin_array[i])
            eu.training_stdev_all[i] = np.std(mainin_array[i])

            #-print normalization values
            print('statistics of ',eu.allfield_names[i],' (training data): mean=',eu.training_mean_all[i],', stdev=',eu.training_stdev_all[i])

    # - normalize input fields wrt spatial mean and stdev (always from training domain, for consistent transformation!)
    for i in range(eu.ninfields):
        mainin_array[i] = (mainin_array[i] - eu.training_mean_all[i])/ eu.training_stdev_all[i]

    mainin_array = ma.masked_invalid(mainin_array)
    mainin_array.filled(0.)
    mainin_array[np.isnan(mainin_array)] = 0.
    print("NAN? ",np.count_nonzero(np.isnan(mainin_array)))

    # - extract selected input fields into single input array
    iuse = 0
    for ifield in range(eu.ninfields):
        if ifield < eu.nmainin:
            bionr_tmp = eu.mainin_fieldnr[ifield]
            # - define used input field arrays
            if bionr_tmp in cf.input_var_use:
                #-array of used bioclim numbers
                eu.mainin_fieldnr_use[iuse] = bionr_tmp
                #-indices of used fields in all fields
                eu.idx_use_in_all[iuse] = ifield
                iuse += 1

        elif ifield == isoilarray:
            #-indices of used fields in all fields
            eu.idx_use_in_all[iuse] = ifield
            iuse += 1

    # - initialize dimensions and arrays for available used input fields
    if use_type == 't':
        eu.ninfields_use = iuse
        eu.ini_global_use()
        print(' eu.ninfields_use=',eu.ninfields_use)

    # - fill array with used fields
    if use_type == 't': print('normalization values of used fields (training data): [name] [mean] [std]')
    data_set = ma.zeros([eu.ninfields_useext, dlat_set, dlon_set])

    for iuse in range(eu.ninfields_use):
        ifield = int(eu.idx_use_in_all[iuse])
        data_set[iuse] = mainin_array[ifield]
        if use_type == 't':
            eu.trainfield_names[iuse] = eu.allfield_names[ifield]
            print(eu.trainfield_names[iuse],eu.training_mean_all[ifield],eu.training_stdev_all[ifield])
            eu.training_mean_use[iuse] = eu.training_mean_all[ifield]
            eu.training_stdev_use[iuse] = eu.training_stdev_all[ifield]

    #-if simpleFit and infields_ext_mode==1: calc and save stats of cross-terms of used fields
    if use_type == 't' and cf.model_training == 'simpleFit' and cf.infields_ext_mode == 1:
        eu.ninfields_useorg = eu.ninfields_use #keep original number of unsed input fields
        icross = eu.ninfields_use
        for juse in range(eu.ninfields_use):
            for iuse in range(juse):
                #-calc cross-terms
                crossfield_tmp = data_set[iuse]*data_set[juse] #from normalized fields(caution: crossfield not normalized!)
                #-save stats and names
                eu.training_mean_use[icross] = np.mean(crossfield_tmp)
                eu.training_stdev_use[icross] = np.std(crossfield_tmp)
                eu.trainfield_names[icross] = eu.trainfield_names[iuse]+' '+eu.trainfield_names[juse]
                icross += 1

    #-if simpleFit and infields_ext_mode==1: normalize cross-terms (for training&investigation data)
    if cf.model_training == 'simpleFit' and cf.infields_ext_mode == 1:
        icross = eu.ninfields_use
        for juse in range(eu.ninfields_use):
            for iuse in range(juse):
                #-calc cross-terms
                crossfield_tmp = data_set[iuse]*data_set[juse] #from normalized fields(caution: crossfield not normalized!)
                #-normalize
                data_set[icross] = (crossfield_tmp - eu.training_mean_use[icross])/ eu.training_stdev_use[icross]
                icross += 1

        #-from now on: treat extended cross-fields as part of infields used
        eu.ninfields_use = eu.ninfields_useext
        print('include cross-terms in used input fields... ->extend ninfields_use=',eu.ninfields_use)

    if use_type == 't': print('indices of used input fields in all input fields(starting with 0):',eu.idx_use_in_all[:])
    data_set.filled(0.)

    iv_set = eu.transpose_data_set(data_set, dlat_set, dlon_set)

    if cf.sample_factor <= 1. or int(cf.sample_factor) != cf.sample_factor:
        if cf.sample_factor != 1:
            print('!!!CAUTION: sample_factor must be a positive intgeger value! Skip sampling...')

        return lat_set, lon_set, dlat_set_out, dlon_set_out, data_set, mainin_array_full, iv_set, mask

    else:

        ### SPATIALLY AVERAGED SUBSAMPLING OF INPUT FIELDS ###
        # Reducing the size of the calculation domain by downsampling to see if that makes a difference
        # Remember to set the radius of the sampling to the correct size equivalently to the downsampling factor!

        # - set reduced dimensions
        sample_dlat_set = dlat_set // cf.sample_factor
        sample_dlon_set = dlon_set // cf.sample_factor
        sample_lat_set = np.zeros(sample_dlat_set)
        sample_lon_set = np.zeros(sample_dlon_set)

        # - get reduced averaged lat-lon and input fields
        sample_mainin_array_full = ma.zeros([eu.ninfields, sample_dlat_set, sample_dlon_set])
        sample_data_set = ma.zeros([eu.ninfields_use, sample_dlat_set, sample_dlon_set])
        for i in range(sample_dlat_set):
            sample_lat_set[i] = np.sum(lat_set[cf.sample_factor*i:cf.sample_factor*i+(cf.sample_factor)])/float(cf.sample_factor)
        for i in range(sample_dlon_set):
            sample_lon_set[i] = np.sum(lon_set[cf.sample_factor*i:cf.sample_factor*i+(cf.sample_factor)])/float(cf.sample_factor)
        for i in range(eu.ninfields):
            for j in range(sample_dlat_set):
                for k in range(sample_dlon_set):
                    sample_mainin_array_full[i,j,k] = ma.mean(mainin_array_full[i,cf.sample_factor*j:cf.sample_factor*j+(cf.sample_factor),cf.sample_factor*k:cf.sample_factor*k+(cf.sample_factor)])
        for i in range(eu.ninfields_use):
            for j in range(sample_dlat_set):
                for k in range(sample_dlon_set):
                    sample_data_set[i,j,k] = ma.mean(data_set[i,cf.sample_factor*j:cf.sample_factor*j+(cf.sample_factor),cf.sample_factor*k:cf.sample_factor*k+(cf.sample_factor)])

        print("sampled latitude points: ",sample_dlat_set)
        print("sampled longitude points: ",sample_dlon_set)

        sample_iv_set = eu.transpose_data_set(sample_data_set, sample_dlat_set, sample_dlon_set)

        return sample_lat_set, sample_lon_set, sample_dlat_set, sample_dlon_set, sample_data_set, sample_mainin_array_full, sample_iv_set, mask


#########################################
##### READ_SITES ########################
# getting archeological site data from file
# - author: Christian Wegener, IGMK, UniKoeln
# - version: ???(from /data/hescor/owf/hep on 2024-12-04)
# - modification history:
#2025-01-21     Annika Vogel    generalize to list of input files
###
def read_sites():
    import math

    ### READ SITE CORRDINATES FROM FILEs ###
    import pandas as pd

    lon_s = []
    lat_s = []
    for path in cf.sites_path:
        # - read table, get lon&lat
        df = pd.read_excel(path)
        print('Reading site data from file:',path)
        lon_file = list(df[cf.sites_lonname])
        lat_file = list(df[cf.sites_latname])

        # - remove invalid entries
        c = 0
        while c < len(lon_file):
            if not isinstance(lon_file[c], float):
                del lon_file[c]
                continue
            if math.isnan(lon_file[c]):
                del lon_file[c]
                continue

            c += 1
        c = 0
        while c < len(lat_file):
            if not isinstance(lat_file[c], float):
                del lat_file[c]
                continue
            if math.isnan(lat_file[c]):
                del lat_file[c]
                continue

            c += 1

        # - put site data from all files into single array
        for i in range(len(lon_file)):
            lon_s.append(lon_file[i])
            lat_s.append(lat_file[i])

        lon_s = np.array(lon_s)
        lat_s = np.array(lat_s)

    '''
    #CAUTION: old, only for preprocessed NCDF files with core ares defined
    else:
        # - read file, select lon&lat
        input_data = Dataset(cf.sites_path)
        if cf.sites_iso_or_not in ['not','no','n','nein','nicht']:
            lat_s = input_data.groups['all_sites'].variables['lat'][:]
            lon_s = input_data.groups['all_sites'].variables['lon'][:]
        else:
            lat_s = input_data.groups['sites_in_isolines'].variables['lat'][:]
            lon_s = input_data.groups['sites_in_isolines'].variables['lon'][:]
        input_data.close()
    '''

    ### optional selection of sub-region ###
    if cf.sites_region == 'west':
        sites = []
        for i in range(len(lat_s)):
            if lon_s[i] <= 10:
                sites.append([lat_s[i], lon_s[i]])
    elif cf.sites_region == 'east':
        sites = []
        for i in range(len(lat_s)):
            if lon_s[i] > 10:
                sites.append([lat_s[i], lon_s[i]])
    else:
        sites = np.zeros([len(lat_s), 2])
        sites[:, 0] = lat_s
        sites[:, 1] = lon_s

    sites = np.array(sites)
    print("Number of read Sites: " + str(len(sites)))

    return sites


#########################################
##### PREP_OUTFILE ######################
# prepare netCDF output file
# - author: Christian Wegener, IGMK, UniKoeln
# - version: ???(from /data/hescor/owf/hep on 2024-12-04)
# - modification history:
#2025-01-20     Annika Vogel    move from main programm into function
#2025-02-28     Annika Vogel    used global variables
###
def prep_outfile(lat,lon,lat_t,lon_t):

    # - define output file and dimensions
    output_data = Dataset(cf.ehep_outpath, 'w', format='NETCDF4_CLASSIC')
    if cf.gridtype_i == 'curvilinear':
        output_data.createDimension('lat', np.shape(lat)[0])
        output_data.createDimension('lon', np.shape(lat)[1])
        output_data.createDimension('lat_t', np.shape(lat_t)[0])
        output_data.createDimension('lon_t', np.shape(lat_t)[1])
    else:
        output_data.createDimension('lat', len(lat))
        output_data.createDimension('lon', len(lon))
        output_data.createDimension('lat_t', len(lat_t))
        output_data.createDimension('lon_t', len(lon_t))

    output_data.createDimension('runs', eu.runs)
    output_data.createDimension('predictors', eu.ninfields_poly)
    # - define output variables
    if cf.gridtype_i == 'curvilinear':
        lat_var_out = output_data.createVariable('lat', np.float32, ('lat', 'lon',))
        lon_var_out = output_data.createVariable('lon', np.float32, ('lat', 'lon',))
        latt_var_out = output_data.createVariable('lat_t', np.float32, ('lat_t', 'lon_t',))
        lont_var_out = output_data.createVariable('lon_t', np.float32, ('lat_t', 'lon_t',))
    else:
        lat_var_out = output_data.createVariable('lat', np.float32, ('lat',))
        lon_var_out = output_data.createVariable('lon', np.float32, ('lon',))
        latt_var_out = output_data.createVariable('lat_t', np.float32, ('lat_t',))
        lont_var_out = output_data.createVariable('lon_t', np.float32, ('lon_t',))

    lat_var_out.longname = "latitude"
    lat_var_out.standard_name = "latitude"
    lat_var_out.units = "degrees_north"
    lon_var_out.longname = "longitude"
    lon_var_out.standard_name = "longitude"
    lon_var_out.units = "degrees_east"

    latt_var_out.longname = "latitude training"
    latt_var_out.standard_name = "latitude"
    latt_var_out.units = "degrees_north"
    latt_var_out.comment = "latitude range of the training set"
    lont_var_out.longname = "longitude training"
    lont_var_out.standard_name = "longitude"
    lont_var_out.units = "degrees_east"
    lont_var_out.comment = "longitude range of the training set"

    hep_var_out = output_data.createVariable('ehep', np.float32, ('runs', 'lat', 'lon',))
    hep_var_out.longname = "Environmental Human Existence Potential"
    auc_var_out = output_data.createVariable('auc', np.float32, ('runs',))
    auc_var_out.longname = "Area Under Receiver Operating Characteristics Curve"
    bs_var_out = output_data.createVariable('bs', np.float32, ('runs',))
    bs_var_out.longname = "Brier Score"
    coef_var_out = output_data.createVariable('coefs', np.float32, ('runs', 'predictors'))
    coef_var_out.longname = "Coefficients of the second degree Logistic Regression based on chosen predictors"
    pre_abs_full_out = output_data.createVariable('pre_abs_full', np.int8, ('lat_t', 'lon_t',))
    pre_abs_full_out.longname = "Presence and absence points full set"
    pre_abs_var_out = output_data.createVariable('pre_abs', np.int8, ('runs', 'lat_t', 'lon_t',))
    pre_abs_var_out.longname = "Presence and absence points used per run"

    return output_data, lat_var_out, lon_var_out, latt_var_out, lont_var_out, hep_var_out, \
            auc_var_out, bs_var_out, coef_var_out, pre_abs_full_out, pre_abs_var_out



#########################################
##### PLOT_PRESABS ######################
# plot 
# - author: Christian Wegener, IGMK, UniKoeln
# - version: ???(from /data/hescor/owf/hep on 2024-12-04)
# - modification history:
#2025-04-24     Annika Vogel    move from ehep_methods.py:pre_abs_sites
###
def plot_presabs(lat,lon,abs_lat,abs_lon,pres_lat,pres_lon,apriabs_lat,apriabs_lon,sites):
    # - load modules
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

    # - initialize plot and map
    fig = plt.figure(figsize=cf.figsize_ref)
    m = Basemap(projection='stere',
                lat_0=(np.min(lat) + np.max(lat)) / 2.,
                lon_0=(np.min(lon) + np.max(lon)) / 2.,
                llcrnrlon = np.min(lon), urcrnrlon = np.max(lon),
                llcrnrlat = np.min(lat), urcrnrlat = np.max(lat),
                resolution='l')
    m.drawcoastlines(linewidth=.6)
    m.drawcountries(linewidth=.4)
    m.drawparallels(np.arange(-50., 100., 10.), linewidth=.3, labels=[True, False, False, False])
    m.drawmeridians(np.arange(0., 360., 10.), linewidth=.3, labels=[False, False, True, False])

    # - plot presence/absence points
    x_1, y_1 = m(abs_lon, abs_lat)
    plt.scatter(x_1, y_1, s=cf.plot_presabs_markersize, c='DarkOrange', edgecolors='none', alpha=.7,label='pseudo-absence')
    x, y = m(pres_lon, pres_lat)
    plt.scatter(x, y, s=cf.plot_presabs_markersize, c='LightSkyBlue', edgecolors='none', alpha=.7,label='presence')
    x_2,y_2 = m(apriabs_lon,apriabs_lat)
    plt.scatter(x_2, y_2, s=cf.plot_presabs_markersize, c='red',edgecolors='none',alpha=.7,label='a-priori absence')

    # - plot sites
    x, y = m(sites[:, 1], sites[:, 0])
    plt.scatter(x, y, s=0.6*cf.plot_presabs_markersize, marker='X',linewidths=0, c='k', label='sites')

    # - inalize plot
    plt.legend(loc='best',framealpha=1.0)
    if cf.annotate:
        plt.annotate(cf.text_anno,xy=(.02,.95),bbox={'facecolor':'white'},xycoords='axes fraction')
    plt.savefig(cf.plot_presabs_path, bbox_inches='tight')
    plt.close(fig)


#########################################
##### PLOT_HISTOGRAM ####################
# plot histogram of normalized fields
# - author: Annika Vogel, IGMK, UniKoeln
# - version: ???(new)
# - modification history:
#2025-04-17     Annika Vogel    move from ehep_methods.py:pre_abs_sites
#2025-05-27     Annika Vogel    overplot gamma-distribution
###
def plot_histogram(pres_indices,abs_indices,apriabs_indices,lat,lon,mainin_array):
    # - load modules
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import ehep_methods as em

    # - prep
    clr_pres = 'LightSkyBlue'
    clr_abs = 'DarkOrange'
    clr_apri = 'red'
    if cf.plot_hist_norm:
        hist_title = 'histograms of normalized input fields'
    else:
        hist_title = 'histograms of input fields'

    hist_ylabel = 'Nr of occurence'
    hist_outname=cf.output_path_common+'/plot_histogram'

    if cf.plot_hist_fieldsel == 'use':
        hist_nsubplots = eu.ninfields_use
        hist_outname=hist_outname+'-use'
    else: #'all'
        hist_nsubplots = eu.ninfields
        hist_outname=hist_outname+'-all'

    if cf.plot_hist_norm:
        hist_outname=hist_outname+'-norm'
        hist_sharex = True
    else:
        hist_sharex = False

    if cf.plot_hist_log:
        hist_outname=hist_outname+'-log'
        hist_title = hist_title+' (log)'
        hist_ylabel = hist_ylabel+' (log)'

    if cf.model_training == 'simpleFit' and cf.plot_hist_norm:
        hist_outname=hist_outname+'-fit'

    # - get all input fields at presence points
    pres_lat_all, pres_lon_all, pres_pred_all = em.fill_labelarrays(pres_indices,lat,lon,mainin_array)
    abs_lat_all, abs_lon_all, abs_pred_all = em.fill_labelarrays(abs_indices,lat,lon,mainin_array)
    apriabs_lat_all, apriabs_lon_all, apriabs_pred_all = em.fill_labelarrays(apriabs_indices,lat,lon,mainin_array)

    # - automatic quadratic arrangement of subplots
    hist_nrows = int(max(np.floor(np.sqrt(hist_nsubplots)),1))
    hist_ncols = int(np.ceil(hist_nsubplots/hist_nrows))
    fig, axs = plt.subplots(hist_nrows, hist_ncols, sharey=True, sharex=hist_sharex, tight_layout=True, figsize=cf.figsize_ref) #hist_nrows/4.*fig.get_size_inches())

    for iplot in range(hist_nsubplots):

        if cf.plot_hist_fieldsel == 'use':
            iall = int(eu.idx_use_in_all[iplot])
            fieldname_tmp = eu.trainfield_names[iplot]
        else:
            iall = iplot
            fieldname_tmp = eu.allfield_names[iplot]

        # - get array
        column_pres = [row[iall] for row in pres_pred_all]
        column_abs = [row[iall] for row in abs_pred_all]
        column_apriabs = [row[iall] for row in apriabs_pred_all]
        #-normalize wrt whole training domain
        if cf.plot_hist_norm:
            column_pres = (column_pres-eu.training_mean_all[iall])/eu.training_stdev_all[iall]
            column_abs = (column_abs-eu.training_mean_all[iall])/eu.training_stdev_all[iall]
            column_apriabs = (column_apriabs-eu.training_mean_all[iall])/eu.training_stdev_all[iall]

        #-pre-calc data limits
        pres_max = np.nanmax(column_pres[:])
        pres_min = np.nanmin(column_pres[:])
        abs_max = np.nanmax(column_abs[:])
        abs_min = np.nanmin(column_abs[:])
        if len(column_apriabs[:]) > 0:
            apri_max = np.nanmax(column_apriabs[:])
            apri_min = np.nanmin(column_apriabs[:])
            data_max = max(pres_max,abs_max,apri_max)
            data_min = min(pres_min,abs_min,apri_min)
        else:
            data_max = max(pres_max,abs_max)
            data_min = min(pres_min,abs_min)
        
        # - plot histogram
        nbins = 20
        if iplot == 0:
            histy,histx,_ = axs.flat[iplot].hist([column_pres,column_abs,column_apriabs], bins=nbins, color=[clr_pres,clr_abs,clr_apri], \
                    label=['pres','pseudo-abs','apriori-abs'], stacked=True, log=cf.plot_hist_log)
        else:
            histy,histx,_ = axs.flat[iplot].hist([column_pres,column_abs,column_apriabs], bins=nbins, color=[clr_pres,clr_abs,clr_apri], \
                    stacked=True, log=cf.plot_hist_log)
        #axs.flat[iplot].plot([0.,0.],color='orange',marker='*',markersize=5)
        
        #-limit x-axis (for normalized data only)
        if cf.plot_hist_norm:
            if data_max > cf.plot_hist_max:
                print('hist: iplot=',iplot,' max exceeded:',abs_max,pres_max) 
                axs.flat[iplot].set_xlim(right=cf.plot_hist_max)
                if len(column_apriabs[:]) > 0:
                    if apri_max > cf.plot_hist_max: axs.flat[iplot].plot(cf.plot_hist_max-0.2,1.5,color='red',alpha=0.8,marker='*',markersize=6)
                if abs_max > cf.plot_hist_max: axs.flat[iplot].plot(cf.plot_hist_max-0.2,1.5,color='orange',alpha=0.8,marker='*',markersize=5)
                if pres_max > cf.plot_hist_max: axs.flat[iplot].plot(cf.plot_hist_max-0.2,1.5,color='blue',alpha=0.8,marker='*',markersize=3)

            if data_min < -cf.plot_hist_max:
                print('hist: iplot=',iplot,' min exceeded:',abs_min,pres_min)                 
                axs.flat[iplot].set_xlim(left=-cf.plot_hist_max)
                if len(column_apriabs[:]) > 0:
                    if apri_min < -cf.plot_hist_max: axs.flat[iplot].plot(-cf.plot_hist_max+0.2,1.5,color='red',alpha=0.8,marker='*',markersize=6)
                if abs_min < -cf.plot_hist_max: axs.flat[iplot].plot(-cf.plot_hist_max+0.2,1.5,color='orange',alpha=0.8,marker='*',markersize=5)
                if pres_min < -cf.plot_hist_max: axs.flat[iplot].plot(-cf.plot_hist_max+0.2,1.5,color='blue',alpha=0.8,marker='*',markersize=3)

        #-plot stats of training domain
        if cf.plot_hist_norm: #by definition of normalization
            mean_all = 0.
            stdev_all = 1.
        else:
            mean_all = eu.training_mean_all[iall]
            stdev_all = eu.training_stdev_all[iall]

        axs.flat[iplot].axvline(x=mean_all,color='orange',linestyle='--') #label = 'mean')
        axs.flat[iplot].axvline(x=mean_all+stdev_all,color='orange', alpha=0.5, linestyle=':')
        axs.flat[iplot].axvline(x=mean_all-stdev_all,color='orange', alpha=0.5, linestyle=':')

        axs.flat[iplot].set_xlabel(fieldname_tmp)
        if cf.plot_hist_norm:
            axs.flat[iplot].xaxis.set_major_locator(ticker.MultipleLocator(1.))
            fig.supxlabel('normalized fields')

        fig.suptitle(hist_title)
        fig.supylabel(hist_ylabel)

        # - calc stats at presence points
        mean_pres = np.nanmean(column_pres)
        std_pres = np.nanstd(column_pres)
        #print(fieldname_tmp,'@pres: mean=',mean_pres,', std=',std_pres)
        #-plot stats at presence points
        axs.flat[iplot].axvline(x=mean_pres,color='blue',linestyle='--') #label = 'mean')
        axs.flat[iplot].axvline(x=mean_pres+std_pres,color='blue', alpha=0.5, linestyle=':')
        axs.flat[iplot].axvline(x=mean_pres-std_pres,color='blue', alpha=0.5, linestyle=':')

        # - define plotting range
        if cf.plot_hist_norm:
            plot_xmin = -cf.plot_hist_max
            plot_xmax = cf.plot_hist_max
        else:
            plot_xmin = histx.min()
            plot_xmax = histx.max()

        # - plot function fitted to presence data by simpleFit
        if cf.model_training == 'simpleFit' and cf.plot_hist_norm:
            from scipy.stats import norm
            dx_bin = (plot_xmax-plot_xmin)/nbins
            x = np.arange(plot_xmin,plot_xmax,dx_bin/10.)
            #x = np.arange(-cf.plot_hist_max,cf.plot_hist_max,0.1)
            gauss = norm(loc=mean_pres,scale=std_pres).pdf(x)
            fit = gauss*len(column_pres)/np.sqrt(2*np.pi*std_pres**2) #scale to area of hist (=number of points*bin width / area of Gaussian)
            #fit = gauss*len(column_pres)*dx_bin/np.sqrt(2*np.pi*std_pres**2)
            axs.flat[iplot].plot(x,fit,color='blue',alpha=0.5)

            # - overplot Gamma distr
            from scipy.stats import skew, gamma
            skew_pres = skew(column_pres, nan_policy='omit')
            if skew_pres < 0.:
                neg_skew = True
                #-flip x-axis (to ensure positive skweness in calc of gamma-distr)
                skew_pres = -skew_pres
                mean_pres = -mean_pres
                x = -x
            else:
                neg_skew = False

            #-calc parameters of distr
            gamma_alpha_pres = 4./skew_pres**2
            gamma_beta_pres = 0.5*std_pres*skew_pres
            gamma_shift_pres = mean_pres -2.*std_pres/skew_pres
            #gamma_shift_pres = -2.*std_pres/skew_pres
            #-get distribution, scale to area
            gamma_fit = gamma.pdf(x,a=gamma_alpha_pres,loc=gamma_shift_pres,scale=gamma_beta_pres)*len(column_pres)/np.sqrt(2*np.pi*std_pres**2)
            if neg_skew: #re-flip x-axis
                x = -x

            #-plot
            axs.flat[iplot].plot(x,gamma_fit,color='black')
           

        #axs.flat[iplot].set_ylim(bottom=1.)
        if cf.plot_hist_log:
            plot_ymax = histy.max()
        else:
            plot_ymax = 200.

        axs.flat[iplot].set_ylim(1.,plot_ymax)
        axs.flat[iplot].set_xlim(plot_xmin,plot_xmax)

    # - finalize
    hist_outname = hist_outname+'.pdf'
    plt.savefig(hist_outname, bbox_inches='tight')
    print('Saving histogram plot as:',hist_outname)
    plt.close(fig)



#########################################
##### PLOT_DISTINCT #####################
# plot distincitveness (all .vs. pres) of human presence conditions for each input field (normalized)
# - author: Annika Vogel, IGMK, UniKoeln
# - version: ???(new)
# - modification history:
#2025-04-21     Annika Vogel    new
#2025-06-24     Annika Vogel    plot total distinctiveness, adjust width to eu.ninfields_use
###
def plot_distinct(pres_means,pres_stds,distinct):
    # - load modules
    import matplotlib.pyplot as plt

    # - prep
    plot_outname=cf.output_path_common+'/plot_distinct.pdf'
    '''
    if cf.plot_hist_fieldsel == 'use':
        hist_nsubplots = eu.ninfields_use
        hist_outname=hist_outname+'-use'
    else: #'all'
        hist_nsubplots = eu.ninfields
        hist_outname=hist_outname+'-all'
    '''
    clr_mean = 'green'
    clr_std = 'purple'
    clr_distinct = 'gray'
    ref_mean = 0.
    ref_std = 1.
    ref_distinct = -np.ceil(np.max(distinct)+.5) #artifically ref for plotting only!
    nplot = len(pres_means)

    # - create plot and plot reference lines for whole training domain
    #if eu.ninfields_use > 100:
    #    fig = plt.figure(figsize=(1.*cf.figsize_ref[0],0.5*cf.figsize_ref[1]))
    #elif eu.ninfields_use > 100
    #else:
    fig = plt.figure(figsize=(0.5*eu.ninfields_use/20.*cf.figsize_ref[0],0.5*cf.figsize_ref[1]))

    #- mean(all training point)=0 (per def of normalization)
    plt.axhline(y=ref_mean,color=clr_mean,linestyle='--',label='mean(all)')
    #- stdev(all training points)=1 (per def of normliazation)
    plt.axhline(y=ref_std,color=clr_std,linestyle='--',label='stdev(all)')
    #- total distinct: artifically ref=-1 for plotting only
    plt.axhline(y=ref_distinct,color=clr_distinct,linestyle='-',label='distinct(base)')

    # - plot presence
    plt.plot(pres_means,'+',color=clr_mean,label='mean(pres)')
    plt.plot(pres_stds,'+',color=clr_std,label='stdev(pres)')
    plt.plot(distinct+ref_distinct,'*',color=clr_distinct,label='distinct')

    # - plot distnctiveness indication (difference presence-to-all)
    for ifield in range(nplot):
        plt.plot([ifield,ifield],[pres_means[ifield],ref_mean],color=clr_mean)
        plt.plot([ifield,ifield],[pres_stds[ifield],ref_std],color=clr_std)
        plt.plot([ifield,ifield],[distinct[ifield]+ref_distinct,ref_distinct],color=clr_distinct)
    
    # - modify plotting labels (names on infields on x-axis)
    plt.xticks(np.arange(nplot), eu.trainfield_names[0:nplot], rotation=60.)

    # - finalize plot
    plt.legend(loc='center left',bbox_to_anchor=(1., 0.5),framealpha=1.0)
    plt.savefig(plot_outname, bbox_inches='tight')
    plt.close(fig)

    


