"""Central configuration for the HEP workflow.
Before running the model data paths need to be defined in this configuration file. 

It defines the experiment setup used by the rest of the pipeline, including:

- Technical runtime options (for example, process count).
- Input data sources and naming conventions (land-sea mask, climate, vegetation, soil, and site files).
- Region/domain bounds and grid metadata for both training and investigation datasets.
- Model preparation and fitting choices (variable selection, presence/absence sampling, run count, and model hyperparameters).
- Plotting and output destinations.

For convenience and easy switching between data and projects several sections are conditionally configured from key strings (for
example, region name in `expname_common` or dataset type in `input_path_t`), so changing those values can automatically switch
related settings.

"""

### Import Modules ###
import copy

### technical config ###
process_count = 8   # Number of processes, should be smaller than cpu number

### general input & config ###
expname_common = 'southern_Africa' #experiment name: common string in land-sea mask & site files !!!CAUTION: #AV-Africa-specific: 'Africa', 'southern_Africa'
# path to input file with land-sea data (string)
path_land_sea_mask = '/data/hescor/anvogel/input-data/topo-data/landmask_30s_'+expname_common+'.nc'


# - inputfield-related setup used for training and investigation (for simplicity)
nbioclim_def = 19           # default number of bioclim variables, required if all in single field (input_onefield_t=true) (default=19, int)
# different specific parts of input files (list of string)
input_filetime_pathfield = ["bio"+"{0:0=2d}".format(i+1) for i in range(nbioclim_def)]
input_filetime_pathend = '_v1.4.0.nc' # common ending of all specific input files (string)
input_filetime_timename = 'time' # name of field in input file indicating time period (for 'time' in input_filedimtype only, string)
input_filetime_timeid = -138000 #-125000 -138000   # EXPERIMENT-SPECIFIC: indicator for selecton of time to be used from input file
                                        #(for 'time' in input_filedimtype only, type dep on fieldname type)

### training input & config (for calculation of HEP parameters) ###
# - grid
gridtype_t = 'lonlat'   # grid type in training input files:  'lonlat'=common lon/lat given for each row/column
                        #                                       'curvilinear'=irregular lon/lat given for each gridpoint
# lat-lon coordinates of domain boundaries for training data (real, in deg N/E ->use neg values for S/W)
if 'southern' in expname_common:
    lat_min_t = -35
    lat_max_t = -16.5
    lon_min_t = 10
    lon_max_t = 40 #sAfrica:exclude desert-population: 17.5
else: #'Africa'
    lat_min_t = -36
    lat_max_t = 40
    lon_min_t = -20
    lon_max_t = 55
# - main input fields (bioclim,vegetation)
# path to main input file for training (string)
input_path_t = 'PATH/TO/INPUT_DATA/bioclim.nc' #e.g. BioClim Dataset by Krapp2021
if 'paleoVeg' in input_path_t: #paleoVeg vegetation fractions
    input_latname_t = 'y'     # name of lat variable in training files (string)
    input_lonname_t = 'x'     # name of lon variable in training files (string)
    input_onefield_t = True   # flag if all input variables are in one field in input file, false=each variable in separate field (flag)
    input_varnames_t = ['vegetation']    # list of names of input fields in training input files (list of string)
    #input_varnames_t = ['Bio1', 'Bio2', 'Bio3','Bio4', 'Bio5', 'Bio6', 'Bio7', 'Bio8','Bio9', 'Bio10', 'Bio11', 'Bio12','Bio13', 'Bio14', 'Bio15',
    #                  'Bio16', 'Bio17', 'Bio18', 'Bio19'] # have to adjust these to allow a list
    input_filedim_type_t = 'field'  # dimension of input fields that is stored in individual input files (string)
    pre_radius_site = 50 # Radius of the presence around site, CAUTION: ~grid resolution (in km, default: 50)
elif 'Krapp' in input_path_t: #Krapp21 bioclim data specific setup
    input_latname_t = 'latitude'     # name of lat variable in training input files (string)
    input_lonname_t = 'longitude'     # name of lon variable in training input files (string)
    input_onefield_t = True   # flag if all variables (eg bioclim/vegetation) are in one field in input file, false=each variable in separate field (flag)
    input_varnames_t = input_filetime_pathfield   # list of names of input fields in training input files (list of string)
    input_filedim_type_t = 'time'  # dimension of input fields that is stored in individual input files (string)
    pre_radius_site = 50 # Radius of the presence around site, CAUTION: ~grid resolution (in km, default: 50)
else: #eg 'Armstrong' #Armstrong bioclim data specific setup
    input_latname_t = 'y'     # name of lat variable in training input files (string)
    input_lonname_t = 'x'     # name of lon variable in training input files (string)
    input_onefield_t = True   # flag if all input variables are in one field in input file, eg as 3rd dimension, false=each variable in separate field (flag)
    input_varnames_t = ['bio_var']    # list of names of input fields in training input files (list of string)
    #input_varnames_t = ['Bio1', 'Bio2', 'Bio3','Bio4', 'Bio5', 'Bio6', 'Bio7', 'Bio8','Bio9', 'Bio10', 'Bio11', 'Bio12','Bio13', 'Bio14', 'Bio15',
    #                  'Bio16', 'Bio17', 'Bio18', 'Bio19'] # have to adjust these to allow a list
    input_filedim_type_t = 'field'  # dimension of input fields that is stored in individual input files (string)
                                # !if 'time': see also 'input_filetime_*' for specific setup!
                                # - TODO'none': only one field at one time in file
                                # - 'field'(default): all fields in one file, only one time in file
                                # - 'time': diff times in one file, only one field in file
                                #       !requires list of biovar-parts of path for each variable in 'input_filetime_pathfield'
                                #       !requires common pre-biovar-part of path in 'biolim_path_t' & post-biovar-part in 'input_filetime_pathend'
                                #       !requires one-entry list 'input_varnames_t' for common fieldname in each file
                                #           OR list with field-dimension of 'input_filetime_pathfield'
                                #       !requires 'input_filetime_{fieldname}&{timeid}'
                                # - TODO'fieldtime' :all fields for diff times in one file, requires time-specific config
    pre_radius_site = 200 # Radius of the presence around site, CAUTION: ~grid resolution (in km, default: 50)

# - soil
# path to input soil file for training (string)
soil_path_t = '/PATH/TO/SOIL_DATA' #'/data/hescor/cwegener/Central_Europe/Band_Neolithikum/data/LBK_soil_map_raster_EU_interpol.nc'
soil_varname_t = 'Band1'    # name of soil field in training soil file (string)

### investigation input & config (for application of HEP parameters) ###
# - grid
gridtype_i = gridtype_t     # grid type in investigation input files ('curvilinear' / 'lonlat')
# lat-lon coordinates of domain boundaries for investigation data (real, in deg)
lat_min_i = lat_min_t
lat_max_i = lat_max_t
lon_min_i = lon_min_t
lon_max_i = lon_max_t ##sAfrica:exclude desert-population:apply to whole sAfrica: 40
#lat_min_i = -36 #TMP:apply to whole Africa...
#lat_max_i = 40 #TMP
#lon_min_i = -20 #TMP
#lon_max_i = 55 #TMP

# - main input fields
# path to input file for training (string)
input_path_i = input_path_t #'/data/hescor/cwegener/Central_Europe/Band_Neolithikum/data/Bioclim_interpol.nc'
input_latname_i = input_latname_t #'y'  # name of lat variable in investigation input files (string)
input_lonname_i = input_lonname_t #'x'  # name of lon variable in investigation input files (string)
input_filedim_type_i = input_filedim_type_t # dimension of input fields that is stored in individual input files (string)
# list of names of input fields in investigation input files (string)
input_varnames_i = copy.copy(input_varnames_t) #['Bio1', 'Bio2', 'Bio3','Bio4', 'Bio5', 'Bio6', 'Bio7', 'Bio8','Bio9', 'Bio10', 'Bio11', 'Bio12','Bio13', 'Bio14', 'Bio15',
#                  'Bio16', 'Bio17', 'Bio18', 'Bio19']
input_onefield_i = input_onefield_t     # flag if all variables are in one field in input file, false=each variable in separate field (flag)
# - soil
soil_path_i = soil_path_t       # path to input soil file for investigation (string)
soil_varname_i = soil_varname_t # name of soil field in investigation soil file (string)


### site input & config ###
sites_region = 'all'           # option for area subselection (default 'all' / 'east':lon>10deg / 'west':lon<=10deg, string)
sites_latname = 'Latitude'     # name of lat variable in site files (string)
sites_lonname = 'Longitude'    # name of lon variable in site files (string)
# path to input files with archeological site data (list of strings)
sites_path = ['/PATH/TO/SITE_DATA'+expname_common+'.xlsx']


### calculation config ###
# - data use
# bioclim using the Number of Bioclim, starting with 1 (not 0 as usual in python!)
input_var_use = list(range(1,nbioclim_def+1))#[1,8,10,16] #default(LBK):[1,2,12,18]
#input_var_use = [1,2,3,4,5,6,7,10,11,12,13,14,16,17,18] #stdev<1-only
#input_var_use = [1,3,4] #paleoVeg-grouped-77ka:stdev<1-only
soil_use = False #flag, if soil data are additionally used (default False, flag)
# select type of limits for apriori absence points: 0=none, 1=predefined 'bio*_min/max', 2=min/max of any pres conditions (for each infield), 3=min/max of all pres cond
absapri_limits_mode = 2
# define limits of bioclim variables for apriori absence points
bio1_min = -30. #minimum limit for annual mean temp [degC]
bio1_max = 160. #maximum limit for annual mean temp [degC]
bio2_min = 0.  #minimum limit for mean diurnal temp range [degC]
bio2_max = 100. #maximum limit for mean diurnal temp range [degC]

# - data preparation
#[input-dependent:defined above] pre_radius_site = 200 # Radius of the presence around site, CAUTION: ~grid resolution (in km, default: 50)
infields_ext_mode = 0   # simpleFit-only: extend input fields: 0=none, bit1=quadratic cross-terms, ...TODO:gradients... (default False, int)
sample_factor = 1 # Sample factor for the downscaling investigation. When altering it, increase the radius similarly (default: 1, CAUTION: integer-only!)
train_absapri_only = False  #flag, if training with apriori absence only, or both pseudo absence and apriori absence (default False, flag)
                            # CAUTION: all absence points are used if too few apriori-absence points for training
train_absapri_minnum = 20 #minimum number of apriori-absence points if training only with them (only for train_absapri_only = True), otherwise: use also pseudo absence (int)
dataratio_trai = 1. #0.8 #relative amount of data to use for training vs testing, default=0.8 (0.<real<=1.)
                        # NOTE: if =1, runs differ only by random selection of pseudo-abs points (if abs_fraction<1)
abs_fraction = 1. #1./3. #relative amount of pseudo-absence data to use, default=1./3. (0.<real<=1.)
                        # NOTE: if =1, runs differ only by random splitting of trainig/test data (if dataratio_trai<1)
#cut_Africa_from_Europe = False #(commented out)

# - calculation
runs = 50 #20 50 ref:1000          # number of different HEP calculations = ensemble size (int)
                                #(=realizations wrt random splitting of training/test data & random selection of pseudo-abs points)
                                #CAUTION: set to 1 if no ranom sampling (if dataratio_trai = 1. & abs_fraction = 1.)
# Choose statistical model to fit training data (string)
model_training = 'logreg'      # 'logreg':logistic regression / 'rf':random forest / 'simpleFit': simple (Gaussian) fit for each input field
                                #(CAUTION: results differ!, default 'logreg', string)
logreg_lasso = 1            # inverse regularization strength of minimization in logistic regression (default 1.0, int)
logreg_tol = 1e-2           # tolerance limit for convergence of logistic regression (defaut 1e-4, real)
logreg_max_iter = 100 #1000 ref:15000  # maximal number of iterations for fit convergence of logistic regression (default 100, int)

### plot & output config ###
output_path_common = '/PATH/TO/OUTPUT'      # common part of output path for plots and data (string)
# - plots
annotate = False        # flag if annotation text to be plotted (flag)
text_anno = "d)"        # annotation text  in plot (string)
figsize_ref = (10,10)   # reference size of plots (tuple of real)

plot_presabs = True     # flag if presence/absence map to be plotted (flag)
plot_presabs_path = output_path_common+'/plot_pres-abs.pdf'     # path to output presence-absence plot (string)
plot_presabs_markersize = 1.*50 #pre_radius_site #ref:160       # markersize in presence-absence plot (real)

plot_hist = True        # flag if histogram of normalized input fields at pres/abs points to be plotted (flag)
plot_hist_fieldsel = 'all' # define which input fields to plot in histogram: 'all' / 'use' (string)
plot_hist_norm = False  # flag if x-values on histogram should be normalized wrt domain statistics(x-mean/stdev) (flag) !CAUTION: fit only for norm!
plot_hist_log = True   # flag if count (y-axis) of histogram should be logaritmic (suggested for small #pres/#abs ratio, default: False flag)
plot_hist_max = 5.      # maximum for normalized x-values to be plotted in histogram

plot_distinct = True    # flag if distinctiveness (all .vs. pres) of human presence conditions to be plotted (flag)

# - data
ehep_outpath = output_path_common+'/hep-out.nc' # path to main output file (string)

