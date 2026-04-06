"""
A land-sea mask is a binary grid that indicates whether each grid cell is land or sea. 

In this case it is synthetically generated for testing purposes, with a simple pattern. 
It serves to demonstrate how to use a land-sea mask within the HEP model.
The land-sea mask should have the same domain as the site locations created with site_locations.py. 

The data that is used for the mask is from Natural Earth Data. (see README for Details)

"""
# UNDER CONSTRUCTION

import os
import sys
import numpy as np
from netCDF4 import Dataset
import cartopy
import cartopy.io.shapereader as shpreader
from shapely.geometry import shape, Point
from shapely.ops import unary_union
from shapely.prepared import prep

# Domain defaults matching configure.py example CAUTION: currently only land, no sea
LAT_MIN = -25.0
LAT_MAX = -15.0
LON_MIN = 15.0
LON_MAX = 25.0


def create_mask_from_naturalearth(lat, lon, cache_dir=None):
    """
	Create a land/sea mask (1=land, 0=sea) by using Natural Earth 'land'.
    lat: 1D array (ny,), lon: 1D array (nx,). Returns ndarray (ny,nx) dtype int8.
    """

    # Acquire Natural Earth 'land' shapefile 
    shp = shpreader.natural_earth(resolution='110m', category='physical', name='land')
    reader = shpreader.Reader(shp)
    geoms = [shape(feature.geometry) for feature in reader.records()]
    if len(geoms) == 0:
        raise RuntimeError("No land geometries found in Natural Earth data")

    land_union = unary_union(geoms)
    land_prep = prep(land_union)

    ny = lat.size
    nx = lon.size
    mask = np.zeros((ny, nx), dtype=np.int8)
    for j in range(ny):
        lat_j = float(lat[j])
        for i in range(nx):
            lon_i = float(lon[i])
            pt = Point(lon_i, lat_j)
            mask[j, i] = 1 if land_prep.contains(pt) or land_prep.touches(pt) else 0

    return mask


def write_mask_to_netcdf(output_path, lat, lon, mask):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with Dataset(output_path, "w", format="NETCDF4") as ds:
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))

        var_lat = ds.createVariable("lat", "f4", ("lat",))
        var_lon = ds.createVariable("lon", "f4", ("lon",))
        var_mask = ds.createVariable("land_sea_mask", "i1", ("lat", "lon"), zlib=True)

        var_lat.units = "degrees_north"
        var_lon.units = "degrees_east"
        var_mask.long_name = "land (1) / sea (0) mask"
        var_mask.flag_values = "0, 1"
        var_mask.flag_meanings = "sea land"

        var_lat[:] = lat
        var_lon[:] = lon
        var_mask[:, :] = mask


if __name__ == "__main__":
    nx = 50
    ny = 50
    lat = np.linspace(LAT_MIN, LAT_MAX, ny)
    lon = np.linspace(LON_MIN, LON_MAX, nx)

    # use cartopy cache inside the project working dir
    cache_dir = os.path.join(os.getcwd(), "cartopy_cache")
    try:
        mask = create_mask_from_naturalearth(lat, lon, cache_dir=cache_dir)
    except Exception as e:
        print("ERROR creating Natural Earth land-sea mask:", str(e), file=sys.stderr)
        sys.exit(2)

    # write results to output directory
    output_path = "/data/hescor/akoepke/HEP-WHB/land_sea_mask.nc"
    write_mask_to_netcdf(output_path, lat, lon, mask)
    print(f"Mask is saved in: {output_path}")







