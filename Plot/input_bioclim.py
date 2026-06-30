"""
This is a plot routine for visualizing bioclimatic variables.
It uses bioclim variables from Krapp et al. (2021).
To use this routine, you need to provide the appropriate input data files 
and specify the desired bioclimatic variable to plot. Furthermore you need to
define the output directory and file name for the plot.

"""



from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import math
from mpl_toolkits.basemap import Basemap

def plot_bioclim():

    # - model input
    bioclim_var = 'bio15'  #bioclim variable to plot, e.g. 'bio01'
    input_pattern = '/PATH/TO/INPUT/DATA/Krapp2021_{var}_v1.4.0.nc'  #path to input variables, '{var}' is replaced by bioclim_var
    lonname = 'longitude'  
    latname = 'latitude'   
    timename = 'time'      #name of time variable in the netCDF file
    time_bp = -51000       #years before present to select from the time dimension
    input_file = input_pattern.format(var=bioclim_var)
    print('reading input file:', input_file)

    # - domain
    lat_min, lat_max = -35.0, -16.5
    lon_min, lon_max = 10.0, 40.0

    # - site input
    plot_presence = True  # flag: overlay presence points from path_sites
    path_sites = ['/PATH/TO/SITES/FILE.xlsx']
    sites_latname = 'Latitude'
    sites_lonname = 'Longitude'

    # - site setup
    site_alpha = 1.0
    site_label = 'presence points'
    site_plotsize = 30
    ignore_sites = not plot_presence

    # - plot config
    colorscheme = None  # None = auto pick temperature/precipitation cmap
    plot_coord_step_lat = 5
    plot_coord_step_lon = 5

    output_dir = '/PATH/TO/OUTPUT/DIRECTORY'
    output_file = output_dir + '/input_' + bioclim_var + '.pdf'

    #####################################

    ds = Dataset(input_file)
    lon = np.array(ds.variables[lonname][:])
    lat = np.array(ds.variables[latname][:])
    var = ds.variables[bioclim_var]
    if timename in var.dimensions and time_bp is not None:
        time = ds.variables[timename][:]
        tidx = int(np.argmin(np.abs(time - time_bp)))
        data = ma.array(var[tidx, :, :])
        actual_time = float(time[tidx])
    else:
        data = ma.array(var[:])
        actual_time = None
    long_name = getattr(var, 'long_name', bioclim_var)
    units = getattr(var, 'units', '')
    ds.close()

    # clip to domain
    lmask = (lon >= lon_min) & (lon <= lon_max)
    amask = (lat >= lat_min) & (lat <= lat_max)
    lon = lon[lmask]
    lat = lat[amask]
    data = data[np.ix_(amask, lmask)]

    if colorscheme is None:
        digits = ''.join(c for c in bioclim_var if c.isdigit())
        colorscheme = 'YlGnBu' if digits and int(digits) >= 12 else 'RdYlBu_r'

    lat_0 = (lat_min + lat_max) / 2.
    lon_0 = (lon_min + lon_max) / 2.

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_axes([0.0, 0.1, 1, 0.8])
    m = Basemap(projection='stere',
                lat_0=lat_0,
                lon_0=lon_0,
                llcrnrlon=lon_min, urcrnrlon=lon_max,
                llcrnrlat=lat_min, urcrnrlat=lat_max,
                resolution='l')
    par = m.drawparallels(np.arange(-90., 90., plot_coord_step_lat), linewidth=.6, labels=[True, False, False, False])
    merid = m.drawmeridians(np.arange(0., 360., plot_coord_step_lon), linewidth=.6, labels=[False, False, True, False])

    lon_2d, lat_2d = np.meshgrid(lon, lat)  # convert 1D lat/lon arrays into 2D arrays
    lon_map, lat_map = m(lon_2d, lat_2d)  # compute map proj coordinates
    pcm = m.pcolormesh(lon_map, lat_map, data, cmap=colorscheme)
    pcm.set_edgecolor('face')

    title = long_name.capitalize() if long_name else bioclim_var
    if actual_time is not None:
        title += ' (' + str(actual_time / 1000) + ' ka BP)'
    fig.suptitle(title, fontsize=14)

    cbar = fig.colorbar(pcm, orientation='horizontal', fraction=0.046, pad=0.04, shrink=0.85)
    cbar.set_label(units if units else bioclim_var)

    if not ignore_sites:
        for path in path_sites:
            print('reading site coordinates from', path)
            df = pd.read_excel(path)
            lon_file = list(df[sites_lonname])
            lat_file = list(df[sites_latname])
            lon_s = []
            lat_s = []
            for lo, la in zip(lon_file, lat_file):
                if not isinstance(lo, float) or math.isnan(lo):
                    continue
                if not isinstance(la, float) or math.isnan(la):
                    continue
                lon_s.append(lo)
                lat_s.append(la)
            sites_lon = np.array(lon_s)
            sites_lat = np.array(lat_s)
            sites_lon_map, sites_lat_map = m(sites_lon, sites_lat)  # compute map proj coordinates
            plt.scatter(sites_lon_map, sites_lat_map, s=site_plotsize, marker='x', c='k',
                        label=site_label, alpha=site_alpha)
        plt.legend(loc='upper right', framealpha=1.0)

    print('=> output plot-file:', output_file)
    plt.savefig(output_file, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    plot_bioclim()
