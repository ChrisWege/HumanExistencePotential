`configure.py`
- configuration of the HEP code

`ehep_run.py`
- main calculation of the HEP

`ehep_inout.py`
- input and output related functions
- reads netcdf file
- reads archaeological site data

`ehep_methods.py`
- fills label-specific index, coordinate, predictor arrays
- specifies presence/pseudo-absense/apriori-absence points
- assigns presence/pseudo-/apiori-absence label to single point
- prepares training and test data
- trains logistic regression, random forest or gaussian fit to data
- does main calculation of HEP (parallel)

`ehep_util.py`
- other utility functions used in the HEP-Model

`plot_routine.py`
- routine for plotting Human Existence Potential

`site_locations.py`
- generates sites locations
- only presence points 

`land_sea_mask.py`
- generates a land-sea-mask using Natural Earth Data


