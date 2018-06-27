##%%writefile jupyter_notebook_exporter_python_mycfg.tpl

# coding: utf-8

#%%

get_ipython().run_line_magic('pylab', 'inline')
# import pickle as pkl
import pandas as pd
from matplotlib import dates
import xarray as xr
# mpl.rcParams['figure.dpi'] = 200

from IPython.display import clear_output
from IPython.core.debugger import Tracer # Tracer()()


#%%

bin_div = pd.read_hdf('pythondata/hiwcdata.h5',key='bin_div').as_matrix().ravel()
bin_mid = (bin_div[:-1]+bin_div[1:])/2.
bin_diff = np.diff(bin_div)
numofflts = 23
rootpath = '/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin'
cmbds = xr.open_dataset('pythondata/rasta_raw_cmb.h5')
rastads = xr.open_dataset('pythondata/rastacombine.h5')
rawpsd = xr.open_dataset('pythondata/psdds.h5')


#%%

filepath = '/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/data_src/data20161012/BOMradar/test'
filename = '20140203_174910.nc'
fullpath = filepath+'/'+filename

ds = xr.open_dataset(fullpath,mask_and_scale=True)
dbz = ds.DBZ
a=np.unique(ds.Coverage.values)
a[~isnan(a)]
dssub = ds.isel(z0=12)
plt.figure()
fig = dssub.DBZ.squeeze().plot.contourf()


#%%

fn = cmbds.flightnum.values


#%%

np.where(np.diff(fn))


#%%

tmpstartind = list(np.where(np.diff(fn))[0])+[-1]


#%%

tmpendind = [0]+list(1+np.where(np.diff(fn))[0])


#%%

str(cmbds.time[tmpendind].values[0])


#%%

len(cmbds.time[tmpendind].values)


#%%

[print('flight '+format(i,'02')+', ',str(x),str(y)) for i,x,y in zip(range(1,24),cmbds.time[tmpendind].values,cmbds.time[tmpstartind].values)]


#%%

cmbds.time[tmpendind].values


#%%

filepath = '/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/data_src/data20161012/BOMradar/test'
filename = '20140203_174910.nc'
fullpath = filepath+'/'+filename

