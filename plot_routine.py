from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import matplotlib.ticker as mticker
import numpy as np
import numpy.ma as ma
import pandas as pd
import math
import cartopy as cpy
import cartopy.feature as cfeat
import cartopy.crs as ccrs
from scipy.interpolate import griddata
import matplotlib.colors as mcolors
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.basemap import Basemap, cm
import itertools

"""Plotting routine for HEP. 

Plots the Human Existence Potential depending on Input Data. 
"""

def plot_potential():

    ############## CONFIG ###############
    # - model input
    stat_type = 'mean' #'mean' 'stdev'
    path = '/data/hescor/akoepke/HEP-WHB'
    expname_common = 'southern_Africa' #'sAfrica' #'southern_Africa'
    hep_indata = 'Krapp21' #'paleoVeg-grouped' #'Krapp21'
    hep_time = '' #'-77ka-E17p5' #'-125ka' #'' '-125ka'
    hep_method = 'logreg' #'simpleFit' 'logreg'
    hep_sampling = 'absFrac1-10' #'50x08' 'det' 'absFrac2-3' 'det_wgt-distinct' 'meanPnegstdev2-det'
    hep_radius = '50' #50 30 80
    hep_infields = 'bioAll' #'bioAll' 'bioAll-excl8-9' 'bio1-8-10-16' 'vegAll'
    hep_other = '_apriLimit2' #'' '_apriLimit2' '_apriLimit2-10' '_apriLimit2-20'
    input_file = path+'/hep-out.nc'
    print('reading intput file:',input_file)

    # - site input
    path_sites = ['/data/hescor/akoepke/HEP-WHB/presence_locations.xlsx'] 
    sites_latname = 'Latitude'     # name of lat variable in site files (string)
    sites_lonname = 'Longitude'

    # - site setup
    site_alpha = 1. #0.3 #0.3 0.7
    site_label = "sites" #old LBK #full LBK
    site_plotsize = 30 #6
    #plot_sites = 'all'    #choose between ['phase1', 'phase2', 'all']
    ignore_sites = False

    # - output config
    stat_longname = 'HEP'
    output_file = path+'/hep-plot.pdf'

    # - plot config
    cbar_location = 'bottom'
    annotate = False
    text_anno = "b)"
    contourbar = True
    plot_title = 'HEP Model idealized'  # e.g. 'My HEP experiment name'
    if stat_type == 'stdev':
        colorscheme = 'YlOrBr' #'YlGnBu' #'viridis' #'RdBu_r' 'seismic' 'bwr' #'coolwarm'
    else:
        colorscheme = 'YlGnBu'

    plot_coord_step_lat = 5
    plot_coord_step_lon = 5
    plot_sea_buffer_boundary = 1 #width of buffer zone in deg which will be additionally plotted at domain boundaries over sea



    #####################################
    #####################################
    #####################################

    input_data = xr.open_dataset(input_file)
    lon = np.array(input_data.lon)
    lat = np.array(input_data.lat)
    accehep = np.ma.array(input_data.variables['ehep'][:][:][:]) #'ehep' #'Acc_HEP'
    
    input_data.close()

    accehep_avg = np.ma.mean(accehep,axis=0)
    accehep_std = np.ma.std(accehep,axis=0)
    
    dlat = len(lat)
    dlon = len(lon)
    lat_flat = lat.flatten()
    lon_flat = lon.flatten()
    points = np.array(list(itertools.product(lat_flat,lon_flat)))
    lat_max = lat.max()
    lat_min = lat.min()
    lon_max = lon.max()
    lon_min = lon.min()
    lat_fine = np.linspace(lat_min,lat_max,401)
    lon_fine = np.linspace(lon_min,lon_max,401)
    grid_lat, grid_lon = np.meshgrid(lat_fine,lon_fine)
    if stat_type == 'stdev' and 'det' not in input_file:
        hep_plot = np.ma.copy(accehep_std)
    else: #'mean'
        hep_plot = np.ma.copy(accehep_avg)
    hep_min = hep_plot.min()
    hep_max = hep_plot.max()
    if stat_type == 'stdev':
        plot_min = 0.
        plot_max = np.ceil(10.*hep_max)/10.
    else:
        plot_min = 0.
        plot_max = 1.

    print('spatial min=',hep_min,', max=',hep_max)
    print('plotting min=',plot_min,', max=',plot_max)
    
    llcrnrlon = lon_min
    urcrnrlon = lon_max
    llcrnrlat = lat_min-plot_sea_buffer_boundary
    urcrnrlat = lat_max
    print('coordinate limits: lat=',lat_min,lat_max,', lon=',lon_min,lon_max)
    if (urcrnrlon-llcrnrlon) > 40:
        plot_coord_step_lon = plot_coord_step_lon*2
    if (urcrnrlat-llcrnrlat) > 20:
        plot_coord_step_lat = plot_coord_step_lat*2
    lat_0 = 50.
    lon_0 = 25.

    fig = plt.figure(figsize=(13, 8.5))
    #print('plot test1')
    ax = fig.add_axes([0.0,0.1,1,0.8]) #, projection=proj)
    m = Basemap(projection = 'stere',
                lat_0 = lat_0,
                lon_0 = lon_0,
                llcrnrlon = llcrnrlon, urcrnrlon = urcrnrlon,
                llcrnrlat = llcrnrlat, urcrnrlat = urcrnrlat,
                resolution = 'l')
    m.drawcoastlines(linewidth=.5, color='Grey')
    m.drawcountries(linewidth=.5, color='Grey')
    par = m.drawparallels(np.arange(-50.,100.,plot_coord_step_lat), linewidth=.6, labels=[True,False,False,False])
    merid = m.drawmeridians(np.arange(0.,360.,plot_coord_step_lon), linewidth=.6, labels=[False,False,True,False])
    lon_2D, lat_2D = np.meshgrid(lon, lat) #convert 1D lat/lon arrays into 2D arrays
    print('hep_plot.shape=',hep_plot.shape)
    lon_map, lat_map = m(lon_2D, lat_2D) # compute map proj coordinates
    POT = m.pcolormesh(lon_map, lat_map, hep_plot, cmap=colorscheme, vmin=plot_min,vmax=plot_max)
    POT.cmap.set_bad('white')
    POT.cmap.set_under('white')
    POT.cmap.set_over('white')
    

    if contourbar == True:
        cbar = fig.colorbar(POT, orientation='horizontal',fraction=0.046,pad=0.04,shrink=0.85, ticks=np.linspace(plot_min,plot_max,11))
        cbar.set_label(stat_longname)

    if annotate == True:
        plt.annotate(text_anno,xy=(.02,.95),bbox={'facecolor':'white'},xycoords='axes fraction')

    POT.set_edgecolor("face")

    sites = 0
    if not ignore_sites:
        for path in path_sites:
            print('reading site coordinates from ',path)
            lon_s = []
            lat_s = []
            df = pd.read_excel(path)
            lon_file = list(df[sites_lonname]) #'x'])
            lat_file = list(df[sites_latname])
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
            for i in range(len(lon_file)):
                lon_s.append(lon_file[i])
                lat_s.append(lat_file[i])
    
            #if sites == 0:
            sites_lon = np.array(lon_s)
            sites_lat = np.array(lat_s)
            '''
            elif sites == 1:
                lon_phase2 = np.array(lon_s)
                lat_phase2 = np.array(lat_s)
                np.delete(lat_phase2, np.where(lon_phase2 > 30))
                np.delete(lon_phase2, np.where(lon_phase2 > 30))
            else:
                lon_core_sites = np.array(lon_s)
                lat_core_sites = np.array(lat_s)
            '''
            sites += 1
    
        sites_lon_map, sites_lat_map = m(sites_lon, sites_lat) # compute map proj coordinates
        plt.scatter(sites_lon_map, sites_lat_map, s=site_plotsize, marker='x', c='k', label=site_label,alpha=site_alpha) #,transform=ccrs.PlateCarree())
    
        plt.legend(loc='upper right',framealpha=1.0)

    print('=> output plot-file:',output_file)
    plt.savefig(output_file, bbox_inches='tight')
    plt.close(fig)

if __name__ == '__main__':
    
    plot_potential()
