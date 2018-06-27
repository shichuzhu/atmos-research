##%%writefile jupyter_notebook_exporter_python_mycfg.tpl

# coding: utf-8

#%%

# %pylab notebook
get_ipython().run_line_magic('pylab', 'inline')
import pandas as pd
import xarray as xr
import glob, re, os
import scipy
# from scipy import io, optimize
import pickle
import datetime
import netCDF4
from mpl_toolkits.basemap import Basemap
# load the notebook import module, refer to $HOME/usr/szhu_setting/python/szpy, set as $PYTHONPATH
from szpy import nbimport
import szpy.sz as sz

## For debug mode
from IPython.core.debugger import Tracer
#Tracer()() #this one triggers the debugger


from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button, RadioButtons
from importlib import reload
# reload(sz)


#%%

# Require 'from szpy import nbimport'
ers = __import__('20170523ellipsoid')
# reload(ers)


#%%

ers.f_plot_E(ers.f_geneigvector(ers.psd))


#%%

# Read all source data file for Darwin campaign
# campaign = 'cayenne'
campaign = 'darwin'
sfsrc, ikpsrc, cdpsrc, robsrc, lampsrc = [ xr.open_dataset('tmp/'+campaign+'/'+campaign+'_sync_bulk.nc',group='/'+x)
                      for x in ['saffire','ikp','cdp','robust','lampproc'] ]

rastanamelist = ['w_ret','Mask_Vz','height_2D','w_wind']
rastasrc = xr.merge( [ xr.open_dataset('tmp/'+campaign+'/'+campaign+'_rasta.nc',group='/'+x)
        for x in rastanamelist ] )

mtsatsrc = xr.open_dataset('tmp/'+campaign+'/'+campaign+'_mtsat.nc',group='/mtsatproc',decode_cf=False)
# Note [] {} are comprehension for lsit and dict, while () is a generator, not a tuple compreh
[ mtsatsrc[x].attrs.pop('missing_value') for x in ['longitude','latitude'] ]
mtsatsrc = xr.conventions.decode_cf(mtsatsrc)

modelsrc = xr.open_dataset('tmp/'+campaign+'/'+'ecmwf.nc')
sdsrc = xr.open_dataset('tmp/'+campaign+'/'+campaign+'_sync_gz.nc',group='/lamp')
bin_div = sdsrc.bin_div
bin_diff = diff(bin_div)
bin_mid = (bin_div[1:] + bin_div[:-1])/2


#%%

# Read all source data file for Cayenne campaign
campaign = 'cayenne'
# campaign = 'darwin'
sfsrc, ikpsrc, cdpsrc, robsrc, lampsrc = [ xr.open_dataset('tmp/'+campaign+'/'+campaign+'_sync_bulk.nc',group='/'+x)
                      for x in ['saffire','ikp','cdp','robust','lampproc'] ]

rastanamelist = ['w_ret','Mask_Vz','height_2D','w_wind']
rastasrc = xr.merge( [ xr.open_dataset('tmp/'+campaign+'/'+campaign+'_rasta.nc',group='/'+x)
        for x in rastanamelist ] )

mtsatsrc = xr.open_dataset('tmp/'+campaign+'/'+campaign+'_mtsat.nc',group='/mtsatproc',decode_cf=True)
# Note [] {} are comprehension for lsit and dict, while () is a generator, not a tuple compreh

modelsrc = xr.open_dataset('tmp/'+campaign+'/'+'ecmwf.nc')
sdsrc = xr.open_dataset('tmp/'+campaign+'/'+campaign+'_sync_gz.nc',group='/lamp')
bin_div = sdsrc.bin_div
bin_diff = diff(bin_div)
bin_mid = (bin_div[1:] + bin_div[:-1])/2


#%%

# sf, ikp, cdp, rob, lamp
get_ipython().system('/sw/netcdf4-4.2-gnu-4.4.6/bin/ncdump -h tmp/darwin/darwin_sync_bulk.nc')


#%%

# rasta
get_ipython().system('/sw/netcdf4-4.2-gnu-4.4.6/bin/ncdump -h tmp/darwin/darwin_rasta.nc')


#%%

# mtsat
get_ipython().system('/sw/netcdf4-4.2-gnu-4.4.6/bin/ncdump -h tmp/darwin/darwin_mtsat.nc')
# haha = netCDF4.Dataset('tmp/darwin/darwin_mtsat.nc')
# { x :haha.groups['mtsatproc'].variables[x].long_name 
#  for x in haha.groups['mtsatproc'].variables.keys()
#  if 'long_name' in haha.groups['mtsatproc'].variables[x].ncattrs() }
# haha.close()


#%%

# model
get_ipython().system('/sw/netcdf4-4.2-gnu-4.4.6/bin/ncdump -h tmp/darwin/ecmwf.nc')


#%%

# routine to select a subset of all data except ECMWF model
# selind = np.where(sfsrc.flightnum.values==13)[0]
selind = np.where((ikpsrc.XKBZR5s.values>.5) & (sfsrc.air_temperature_rm.values<-15)  )[0]
#                  & (lampsrc.validbinnum.values>300) )[0]
sf, ikp, cdp, rob, lamp, rasta, mtsat, sd = [ x.isel(timeutc=selind)
    for x in [sfsrc, ikpsrc, cdpsrc, robsrc, lampsrc, rastasrc, mtsatsrc, sdsrc] ]


#%%

get_ipython().run_cell_magic('capture', '', "%matplotlib notebook\n# Required for the next cell\nwith open('tmp/darwinphasefig.p','rb') as f:\n    axb = pickle.load(f)[1]")


#%%

# interactive plot of probing the phase space
get_ipython().run_line_magic('matplotlib', 'notebook')
# %matplotlib inline

ax = plt.subplot(111, projection='3d')
plt.subplots_adjust(left=0.25, bottom=0.3)
fig = plt.gcf()
x,y,z = lamp.nml.values.T
ax.scatter(x,y,z)
ax.set_xlabel('N0')
ax.set_ylabel('mu')
ax.set_zlabel('lambda')
# The importance of making a copy here. Otherwise xl changes as the ax is changed.
# xl = ax.get_xlim().copy()
# yl = ax.get_ylim().copy()
# zl = ax.get_zlim().copy()
xl = axb.get_xlim().copy()
yl = axb.get_ylim().copy()
zl = axb.get_zlim().copy()

axcolor = 'lightgoldenrodyellow'
axn0 = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axmu = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
axld = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)

sn0 = Slider(axn0, 'n0', -100., 100., valinit=0.)
smu = Slider(axmu, 'mu', -100., 100., valinit=0.)
sld = Slider(axld, 'lambda', -100., 100., valinit=0.)

def update(wtfisthis):
    global sn0, smu, sld
    ax.set_xlim(xl+(xl[1]-xl[0])*sn0.val/100)
    ax.set_ylim(yl+(yl[1]-yl[0])*smu.val/100)
    ax.set_zlim(zl+(zl[1]-zl[0])*sld.val/100)
    fig.canvas.draw_idle()
sn0.on_changed(update)
smu.on_changed(update)
sld.on_changed(update)
update(None)

def on_key(event):
    global xl,yl,zl
    if event.key == 'right':
        # The importance of making a copy here. Otherwise xl changes as the ax is changed.
        xl = ax.get_xlim().copy()
        yl = ax.get_ylim().copy()
        zl = ax.get_zlim().copy()
        [ x.set_val(0) for x in [sn0,smu,sld] ]
    elif event.key == 'left':
        [ x.set_val(0) for x in [sn0,smu,sld] ]
cid = fig.canvas.mpl_connect('key_press_event', on_key)

plt.show()


#%%

sf.timeutc


#%%

ikp


#%%

condition = (sf.flightnum==13) & (sf.timeutc>np.datetime64('2014-02-03T05:43')) & (ikp.XKBZR5s>2.2) & (ikp.Slong<129.75) &(ikp.Slong>129.55) & (sf.upward_air_velocity>0)
ind = sf.timeutc[condition]

tmp = sf.sel(timeutc=ind)
tmpikp = ikp.sel(timeutc=ind)
tmpsd = sd.sel(timeutc=ind)

psd = tmpsd.psddmax.values

get_ipython().run_line_magic('matplotlib', 'inline')
for hehe in psd:
    plt.loglog(bin_mid,hehe,'r')

condition = (sf.flightnum==13) & (sf.timeutc>np.datetime64('2014-02-03T05:43')) & (ikp.XKBZR5s>2.2) & (ikp.Slong<129.75) &(ikp.Slong>129.55) & (sf.upward_air_velocity<0)
ind = sf.timeutc[condition]

tmp = sf.sel(timeutc=ind)
tmpikp = ikp.sel(timeutc=ind)
tmpsd = sd.sel(timeutc=ind)

psd = tmpsd.psddmax.values

for hehe in psd:
    plt.loglog(bin_mid,hehe,'b',alpha=.5)
# Red lines are in the updraft core, while blue lines are in the downdraft region

