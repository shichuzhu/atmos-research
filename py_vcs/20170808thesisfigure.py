##%%writefile jupyter_notebook_exporter_python_mycfg.tpl

# coding: utf-8

#%%

# %pylab notebook
get_ipython().run_line_magic('pylab', 'inline')
import pandas as pd
import xarray as xr
import glob, re, os
import scipy
import pickle
import datetime
import netCDF4
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
# load the notebook import module, refer to $HOME/usr/szhu_setting/python/szpy, set as $PYTHONPATH
from szpy import nbimport
import szpy.sz as sz

## For debug mode
from IPython.core.debugger import Tracer
#Tracer()() #this one triggers the debugger

# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.widgets import Slider, Button, RadioButtons
# from importlib import reload


#%%

# Change the working directory to hiwc folder
workdir = os.getcwd()
print(os.getcwd())
# For windows home machine.
os.chdir('D:\shared\gmaildrive\hiwc')
print(os.getcwd())

# Require 'from szpy import nbimport'
ers = __import__('20170523ellipsoid')
# reload(ers)


#%%

ers.f_plot_E(ers.f_geneigvector(ers.psd))


#%%

# Read all source data file for Darwin campaign
# campaign = 'cayenne'
campaign = 'darwin'
sfsrc, ikpsrc, cdpsrc, robsrc, lampsrc = [
    xr.open_dataset(
        'tmp/' + campaign + '/' + campaign + '_sync_bulk.nc', group='/' + x)
    for x in ['saffire', 'ikp', 'cdp', 'robust', 'lampproc']
]

rastanamelist = ['w_ret', 'Mask_Vz', 'height_2D', 'w_wind']
rastasrc = xr.merge([
    xr.open_dataset(
        'tmp/' + campaign + '/' + campaign + '_rasta.nc', group='/' + x)
    for x in rastanamelist
])

mtsatsrc = xr.open_dataset(
    'tmp/' + campaign + '/' + campaign + '_mtsat.nc',
    group='/mtsatproc',
    decode_cf=False)
# Note [] {} are comprehension for lsit and dict, while () is a generator, not a tuple compreh
[mtsatsrc[x].attrs.pop('missing_value') for x in ['longitude', 'latitude']]
mtsatsrc = xr.conventions.decode_cf(mtsatsrc)
mtsatsrc['index'] = xr.DataArray(
    np.array(range(mtsatsrc.dims['timeutc'])), dims=['timeutc'])

modelsrc = xr.open_dataset('tmp/' + campaign + '/' + 'ecmwf.nc')
sdsrc = xr.open_dataset(
    'tmp/' + campaign + '/' + campaign + '_sync_gz.nc', group='/lamp')
bin_div = sdsrc.bin_div
bin_diff = diff(bin_div)
bin_mid = (bin_div[1:] + bin_div[:-1]) / 2


#%%

# Read all source data file for Cayenne campaign
campaign = 'cayenne'
# campaign = 'darwin'
sfsrc, ikpsrc, cdpsrc, robsrc, lampsrc = [
    xr.open_dataset(
        'tmp/' + campaign + '/' + campaign + '_sync_bulk.nc', group='/' + x)
    for x in ['saffire', 'ikp', 'cdp', 'robust', 'lampproc']
]

rastanamelist = ['w_ret', 'Mask_Vz', 'height_2D', 'w_wind']
rastasrc = xr.merge([
    xr.open_dataset(
        'tmp/' + campaign + '/' + campaign + '_rasta.nc', group='/' + x)
    for x in rastanamelist
])

mtsatsrc = xr.open_dataset(
    'tmp/' + campaign + '/' + campaign + '_mtsat.nc',
    group='/mtsatproc',
    decode_cf=True)
# Note [] {} are comprehension for lsit and dict, while () is a generator, not a tuple compreh
mtsatsrc['index'] = xr.DataArray(
    np.array(range(mtsatsrc.dims['timeutc'])), dims=['timeutc'])

modelsrc = xr.open_dataset('tmp/' + campaign + '/' + 'ecmwf.nc')
sdsrc = xr.open_dataset(
    'tmp/' + campaign + '/' + campaign + '_sync_gz.nc', group='/lamp')
bin_div = sdsrc.bin_div
bin_diff = diff(bin_div)
bin_mid = (bin_div[1:] + bin_div[:-1]) / 2


#%%

# sf, ikp, cdp, rob, lamp
get_ipython().system('ncdump -h tmp/darwin/darwin_sync_bulk.nc')


#%%

# rasta
get_ipython().system('ncdump -h tmp/darwin/darwin_rasta.nc')


#%%

# mtsat
get_ipython().system('ncdump -h tmp/darwin/darwin_mtsat.nc')


#%%

# model
get_ipython().system('ncdump -h tmp/darwin/ecmwf.nc')


#%%

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'purple', 'brown']
linestyle = ['--', '-']
kws = [{'color': x, 'linestyle': y} for y in linestyle for x in colors]

lonlat = [100, 150, -25, 0]
# plt.figure(figsize=(11, 11))
m = Basemap(
    llcrnrlon=lonlat[0],
    llcrnrlat=lonlat[2],
    urcrnrlon=lonlat[1],
    urcrnrlat=lonlat[3],
    projection='merc',
    fix_aspect=True,
    resolution=None)
m.bluemarble()
m.drawmeridians(np.arange(lonlat[0], lonlat[1] + 0.1, 10), labels=[1, 1, 1, 1])
m.drawparallels(np.arange(lonlat[2], lonlat[3] + 0.1, 5), labels=[1, 1, 1, 1])

mlon, mlat = m(sfsrc.longitude.values, sfsrc.latitude.values)
for szi in range(1, 24):
    tmpind = (sfsrc.flightnum.values == szi)
    tmpiwc = ikpsrc.XKBZR5s.values[tmpind] > 1.5
    if not any(tmpiwc):
        print(str(szi) + ' Emtpy!')
        continue
    kw = kws.pop()
    plt.plot(mlon[tmpind], mlat[tmpind], label=str(szi), lw=2, **kw)

plt.legend(ncol=2)
plt.tight_layout()
fig = plt.gcf()

fig.savefig(
    workdir + '/fig/darwin_flt.pdf', transparent=True, bbox_inches='tight')


#%%

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'purple', 'brown']
colors = colors[0:8]
linestyle = ['--', '-']
kws = [{'color': x, 'linestyle': y} for y in linestyle for x in colors]

lonlat = [-70, -40, 0, 15]
# plt.figure(figsize=(11, 11))
m = Basemap(
    llcrnrlon=lonlat[0],
    llcrnrlat=lonlat[2],
    urcrnrlon=lonlat[1],
    urcrnrlat=lonlat[3],
    projection='merc',
    fix_aspect=True,
    resolution=None)
m.bluemarble()
m.drawmeridians(np.arange(lonlat[0], lonlat[1] + 0.1, 10), labels=[1, 1, 1, 1])
m.drawparallels(np.arange(lonlat[2], lonlat[3] + 0.1, 5), labels=[1, 1, 1, 1])

mlon, mlat = m(sfsrc.longitude.values, sfsrc.latitude.values)
for szi in range(9, 27):
    tmpind = (sfsrc.flightnum.values == szi)
    tmpiwc = ikpsrc.XKBZR5s.values[tmpind] > 1.5
    if not any(tmpiwc):
        print(str(szi) + ' Emtpy!')
        continue
    kw = kws.pop()
    plt.plot(mlon[tmpind], mlat[tmpind], label=str(szi), lw=2, **kw)

plt.legend(ncol=2)
plt.tight_layout()
fig = plt.gcf()

fig.savefig(
    workdir + '/fig/cayenne_flt.pdf', transparent=True, bbox_inches='tight')


#%%

# routine to select a subset of all data except ECMWF model
# selind = np.where(sfsrc.flightnum.values==13)[0]
selind = np.where((ikpsrc.XKBZR5s.values > .5) &
                  (sfsrc.air_temperature_rm.values < -15))[0]
#                  & (lampsrc.validbinnum.values>300) )[0]
sf, ikp, cdp, rob, lamp, rasta, mtsat, sd = [
    x.isel(timeutc=selind)
    for x in
    [sfsrc, ikpsrc, cdpsrc, robsrc, lampsrc, rastasrc, mtsatsrc, sdsrc]
]

