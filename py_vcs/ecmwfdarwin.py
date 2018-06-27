#!/data/keeling/a/szhu28/usr/anaconda2/envs/py36/bin/python
import warnings
# warnings.filterwarnings('ignore')
from numpy import *
import numpy as np
import pandas as pd
import xarray as xr
import glob, re, os
import scipy
from scipy import io, optimize
import pickle
import datetime
import netCDF4
import pygrib
import resource
## For debug mode
from IPython.core.debugger import Tracer
#Tracer()() #this one triggers the debugger

# Convert and create NCAR ECMWF archive
import pygrib
import resource

# Parameters for Darwin
filepaths = glob.glob('/data/gpm/a/shared/szhu28/hiwcproc/era_interim/2014**/*regn128sc*',recursive=True)
kw = {'lat1':-20,'lat2':10,'lon1':110,'lon2':160}
outputnoext = 'tmp/darwin/ecmwf'

filename = '/net/san-b8-ib/data/gpm/a/shared/szhu28/hiwcproc/era_interim/201402/ei.oper.an.pl.regn128sc.2014020306'
# with pygrib.open(filename) as grbs:
grbs = pygrib.open(filename)
grbs.rewind()
varfullnames = {x.shortName:(x.name, x.units) for x in grbs}
numofvar = len(varfullnames)
grbs.rewind()
varnames = [x.shortName for x in grbs[0:numofvar]]

# [ (i,varnames[i],varfullnames[varnames[i]]) for i in range(numofvar) ]

filename = '/net/san-b8-ib/data/gpm/a/shared/szhu28/hiwcproc/era_interim/201402/ei.oper.an.pl.regn128uv.2014020306'
# with pygrib.open(filename) as grbs:
grbsuv = pygrib.open(filename)
grbsuv.rewind()
varfullnamesuv = {x.shortName:(x.name, x.units) for x in grbsuv}
numofvaruv = len(varfullnamesuv)
grbsuv.rewind()
varnamesuv = [x.shortName for x in grbsuv[0:numofvaruv]]

# [ (i,varnamesuv[i],varfullnamesuv[varnamesuv[i]]) for i in range(numofvaruv) ]

grbs.rewind()
level = np.array([x.level for x in grbs[0::numofvar] ]).astype(float)

key = 'level'
dimlvl = xr.DataArray(level, dims=['level'], coords={'level':level}, attrs={
    'units':'hPa'}, name=key)

_dump, lats, lons = grbs[1].data(**kw)
lat, lon = lats[:,0], lons[0,:]
key = 'lat'
dimlat = xr.DataArray(lat, dims=['lat'], coords={'lat':lat}, attrs={
    'units':'degree north'}, name=key)
key = 'lon'
dimlon = xr.DataArray(lon, dims=['lon'], coords={'lon':lon}, attrs={
    'units':'degree east'}, name=key)

encoding = dict(zip([varfullnames[x][0] for x in varnames],[{'zlib':True, 'complevel':1}]*len(varnames)))
encodinguv = dict(zip([varfullnamesuv[x][0] for x in varnamesuv],[{'zlib':True, 'complevel':1}]*len(varnamesuv)))
encoding = {**encoding,**encodinguv}

def f_read_model_pair(filepair):
    grbs = pygrib.open(filepair[0])
    grbsuv = pygrib.open(filepair[1])
    timemod = str(grbs[1].dataDate)+"{:04d}".format(grbs[1].dataTime)
    tmponefile = {}
    for skey in varnames:
        key = varfullnames[skey][0]
        units = varfullnames[skey][1]
        tmpdata = np.array([ x.data(**kw)[0] for x in grbs[varnames.index(skey)::numofvar] ])
        tmponefile[key] = xr.DataArray(tmpdata, coords=[dimlvl,dimlat,dimlon], attrs={'units':units}, name=key)
    for skey in varnamesuv:
        key = varfullnamesuv[skey][0]
        units = varfullnamesuv[skey][1]
        tmpdata = np.array([ x.data(**kw)[0] for x in grbsuv[varnamesuv.index(skey)::numofvaruv] ])
        tmponefile[key] = xr.DataArray(tmpdata, coords=[dimlvl,dimlat,dimlon], attrs={'units':units}, name=key)

    onefile = xr.Dataset(tmponefile, coords={'level':dimlvl,'lat':dimlat,'lon':dimlon}, attrs={'timemod':timemod})
    return onefile, timemod

filepaths.sort()
filepairs = [ ( x, re.sub('regn128sc','regn128uv',x) ) for x in filepaths ]

dslist = []
timestrlist = []
i=0
for x in filepairs:
    print(i, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    i+=1
    ds, timestr = f_read_model_pair(x)
    dslist.append(ds)
    timestrlist.append(timestr)
    
with open(outputnoext+'_tmpstore.p','bw') as f:
    pickle.dump([dslist,timestrlist],f)
    
timemods = pd.DatetimeIndex(timestrlist)
timemod = xr.DataArray(timemods, dims=['timemod'], coords={'timemod':timemods}, name='timemod')
xr.concat(dslist,dim=timemod).to_netcdf(outputnoext+'.nc',mode='w',
    format='NETCDF4',engine='netcdf4',group='/', encoding=encoding)
