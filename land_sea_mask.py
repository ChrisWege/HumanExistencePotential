"""
A land-sea mask is a binary grid that indicates whether each grid cell is land or sea. 
In this case it is synthetically generated for testing purposes, with a simple pattern looking like an island. 
It is not based on real-world data, but serves to demonstrate use a land-sea mask within the HEP model.
The land-sea mask is within a 50x50 grid and can be changed depending on use case. It should have the same domain 
as the site locations created with site_locations.py. 

"""

# UNDER CONSTRUCTION

import numpy as np
from netCDF4 import Dataset


def create_island_mask(nx=50, ny=50, island_width=0.55, island_height=0.65):
	# Create 2D binary mask with an elliptical island centered in the domain
	x = np.linspace(-1.0, 1.0, nx) # 1D array of x-coordinates from -1 to 1
	y = np.linspace(-1.0, 1.0, ny) # 1D array of y-coordinates from -1 to 1
	xx, yy = np.meshgrid(x, y) 

	# Ellipse equation: values inside are land (1), outside are sea (0).
	ellipse = (xx / island_width) ** 2 + (yy / island_height) ** 2 <= 1.0
	return ellipse.astype(np.int8)


def write_mask_to_netcdf(output_path, lat, lon, mask):
	# Write mask to NetCDF with lat/lon dimensions and land_sea_mask variable
	with Dataset(output_path, "w", format="NETCDF4") as ds:
		ds.createDimension("lat", len(lat))
		ds.createDimension("lon", len(lon))

		var_lat = ds.createVariable("lat", "f4", ("lat",))
		var_lon = ds.createVariable("lon", "f4", ("lon",))
		var_mask = ds.createVariable("land_sea_mask", "i1", ("lat", "lon"))

		var_lat.units = "degrees_north"
		var_lon.units = "degrees_east"
		var_mask.long_name = "synthetic land sea mask for testing"
		var_mask.flag_values = "0, 1"
		var_mask.flag_meanings = "sea land"

		var_lat[:] = lat
		var_lon[:] = lon
		var_mask[:, :] = mask




