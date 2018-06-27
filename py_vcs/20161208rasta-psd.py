##%%writefile jupyter_notebook_exporter_python_mycfg.tpl

# coding: utf-8
<script>
  jQuery(document).ready(function($) {  
  $(window).load(function(){
    $('#preloader').fadeOut('slow',function(){$(this).remove();});
  });
  });
</script>
<style type="text/css">
  div#preloader { position: fixed; 
      left: 0; 
      top: 0; 
      z-index: 999; 
      width: 100%; 
      height: 100%; 
      overflow: visible; 
      background: #fff url('http://preloaders.net/preloaders/720/Moving%20line.gif') no-repeat center center; 
  }
</style>
<div id="preloader">
</div>

<script>
  function code_toggle() {
    if (code_shown){
      $('div.input').hide('500');
      $('#toggleButton').val('Show Code')
    } else {
      $('div.input').show('500');
      $('#toggleButton').val('Hide Code')
    }
    code_shown = !code_shown
  } 
  
  $( document ).ready(function(){
    code_shown=false; 
    $('div.input').hide()
  });
</script>
<form action="javascript:code_toggle()"><input type="submit" id="toggleButton" value="Show Code"></form>
#%%

## Cell required to run other cells in the file
## Preferably run under py3k env
get_ipython().run_line_magic('pylab', 'inline')
import pandas as pd
import h5py
import xarray as xr
import glob, re, os
import scipy
import pickle
import sz
from IPython.core.debugger import Tracer

import traceback
import warnings
import sys
# warnings.filterwarnings('error')

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    traceback.print_stack()
    log = file if hasattr(file,'write') else sys.stderr
    log.write(warnings.formatwarning(message, category, filename, lineno, line))
## uncomment the following line to have detailed warning info
# warnings.showwarning = warn_with_traceback

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font',**{'family':'Bitstream Vera Sans','sans-serif':['Helvetica']})

### Load the created data
rawvar = pd.read_hdf('pythondata/hiwcdata.h5',key='rawvar')
rawvards = xr.open_dataset('pythondata/rawvards.h5')
rawpsd = pd.read_hdf('pythondata/hiwcdata.h5',key='rawpsd')
bin_div = pd.read_hdf('pythondata/hiwcdata.h5',key='bin_div').as_matrix().ravel()

# file not yet opened
with open('pythondata/mmd.p', 'rb') as f:
    # file opened
    MMD = pickle.load(f)
    # file opened
# file closed

with open('pythondata/fitted.p', 'rb') as f:
    # file opened
    FITS = pickle.load(f)
    # file opened
# file closed

numofflts = 23
rootpath = '/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin'

## Creating iwc
iwc = rawvar['XKBZR5s']
temp = rawvar['air_temperature']
## Creating NML
targetarr = []
for x in FITS:
    ## None means less than 10 valid bins are present,
    ## 'Error' means the fitting fails for the PSD
    if (x is None) or (x == 'Error'):
        targetarr.append([np.nan]*3)
    else:
        targetarr.append(x.x)
#         print('Warning! Not aligned!')
#         print(x)
NML = np.array(targetarr)

## Create dataset dataf combining rawvar and nml, which is still same len as rawvar
temp = rawvar['air_temperature']
nml = np.empty((len(iwc),3))+np.nan
nml[rawpsd['indpsdbackward'],:] = NML
n0 = nml[:,0]
mu = nml[:,1]
ld = nml[:,2]
mmd = np.empty((len(iwc),))+np.nan
mmd[rawpsd['indpsdbackward']] = MMD
indpsdforward = rawvar['indpsdforward']

dataf = rawvar
dataf['temp'] = temp
dataf['indpsdforward'] = indpsdforward
dataf['n0'] = n0
dataf['mu'] = mu
dataf['ld'] = ld
dataf['mmd'] = mmd
dataf['iwc'] = iwc

### add the psd characteristics to rawvards as well
rawvards['mmd']=xr.DataArray(mmd,coords={'time':rawvards.time})
rawvards['n0']=xr.DataArray(n0,coords={'time':rawvards.time})
rawvards['mu']=xr.DataArray(mu,coords={'time':rawvards.time})
rawvards['ld']=xr.DataArray(ld,coords={'time':rawvards.time})

def loadRastaflt(szi):
    datapath='/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/RASTA/data/'
    rastafn=glob.glob(datapath+'*_F'+str(szi)+'_*.nc')
    if len(rastafn) ==0:
        return None
    ds = xr.open_dataset(rastafn[0])
    tmp = (ds.time.values*3600).astype('timedelta64[s]')
    midnight = np.datetime64(rawvar.index[rawvar['flightnum']==szi][0].date())
    ds['timeSec'] = ds.time
    ds['time'] = xr.DataArray(tmp+midnight,coords={'time':ds.time})
    tmp = ds.time
    if tmp.shape == np.unique(tmp).shape:
        print('flight '+str(szi)+' good')
    else:
        print(str(tmp.shape - np.unique(tmp).shape)+' data duplicate disregarded in the future')
    return ds

##future indexing method: ds.sel(time=datetimearray)

## Easy function to sub-setting data based on iwc, temperature, mmd based on dataf
def subsetdata(tempr=[-inf,inf], iwcr=[-inf,inf], mmdr=[-inf,inf],msd=False):

    tmpind = (dataf['temp']>tempr[0]) & (dataf['temp']<tempr[1])
    tmpind = tmpind & (dataf['iwc']>iwcr[0]) & (dataf['iwc']<iwcr[1])
    # tmpind = (dataf['temp']>-37.5) & (dataf['temp']<-32.5)
    # tmpind = tmpind & (dataf['iwc']<5) & (dataf['iwc']>2.5)
    tmpind = tmpind & (dataf['mmd']>mmdr[0]) & (dataf['mmd']<mmdr[1])
    tmpdataf = dataf[tmpind]
    
    tmpindrawvar = tmpdataf['indpsdforward']
    tmpindrawvar = tmpindrawvar[~isnan(tmpindrawvar)].astype(int)
    
    tmpindrawvar = rawpsd['indpsdbackward'][tmpindrawvar]

    indlvl1 = ~isnan(tmpdataf['indpsdforward'])
    a = tmpdataf['indpsdforward'][indlvl1]
    nonpsd = tmpdataf[indlvl1]
    if msd==False:
        psd = rawpsd['psddmax'].iloc[ a.astype(int),: ]
    else:
        psd = rawpsd['msddmax'].iloc[ a.astype(int),: ]
    ## Use rawvar.iloc[tmpindrawvar] for data access
    return nonpsd,psd, tmpindrawvar

## Merge the rasta dataset with the rawvar dataset
rawvar.index.rename('time',inplace=True)
'''
rastads = xr.open_dataset('rastacombine.h5')
rawds = rawvar.to_xarray()
## Find duplicate conflict variables and rename the rasta ones
rawkeys = set(rawds.data_vars.keys())
rastakeys = set(rastads.data_vars.keys())
dup = rawkeys.intersection(rastakeys)
renamedict = { dupstr:dupstr+'_rasta' for dupstr in dup}
rastads.rename(renamedict,inplace=True)
rastaraw = xr.merge([rawds,rastads])

# rastaraw.to_netcdf('rasta_raw_cmb.h5',format='NETCDF4',mode='w')
'''
## Read rasta rawvar combined xarray dataset
cmbds = xr.open_dataset('pythondata/rasta_raw_cmb.h5')
rastads = xr.open_dataset('pythondata/rastacombine.h5')
rawpsd = xr.open_dataset('pythondata/psdds.h5')


#%%

## Method ONE to integrate different dataset -- mostly recommended
plt.rcParams['image.cmap'] = 'jet' # Set default colormap to jet
mpl.rc("savefig", dpi=400)
plt.rcParams.update({'font.size': 15}) # Set default fontsize to 20 to accomodate figsize=(11,*)
fig=plt.figure(figsize=[11,11])    ########### comment me

for pl,i in zip([1,2,3],[12,13,16]):
    ind = (rawvards.flightnum==i)
    raw = rawvards.isel(time=np.squeeze(np.where(ind)))

    ### plt.contourf won't work if the X or Y contains any NaN!!!!
    ### Therefore need to get rid of all the NaN's in the height_2D array produced by the interpolating.
    ind = all(~isnan(rastads.height_2D),axis=1)
    rasta = rastads.isel(time=np.squeeze(np.where(ind)))

    rawp, rastap = xr.align(raw,rasta,join='inner')

    tmp=rastap['Z']
    # tmp=tmp.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
    # tmp=tmp.where((rastap.Mask_Vz==1))
    tmp=tmp.transpose()

    x = np.tile(rawp.time,(500,1))
    y = rastap.height_2D.transpose().values
    z = tmp
    ax=plt.subplot(3,1,pl)
    plt.contourf(x,y,z,np.linspace(-10,30,11),extend='both')

    time=rawp.time
    alt=rawp.altitude/1e3
    plt.plot(time,alt,'k--',linewidth=1.)

#     ax.autofmt_xdate()
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))

    # ax=plt.gca()
    hcb=plt.colorbar()
    hcb.set_label('Z dBZ')
    plt.ylim(0,20)
    plt.ylabel('Altitude (km)')
    plt.title('F'+str(i))

    ax2=ax.twinx()
    time=rawp.time
    iwc=rawp.XKBZR5s
    ax2.plot(time,iwc,'m',linewidth=1.)
    ax2.set_ylim(0,10)
    _=ax2.set_ylabel('IWC g/m3')
    ax2.yaxis.set_label_coords(.97, .5)


#%%

## Method ONE to integrate different dataset -- mostly recommended
ind = (rawvards.flightnum==12)
raw = rawvards.isel(time=np.squeeze(np.where(ind)))

### plt.contourf won't work if the X or Y contains any NaN!!!!
### Therefore need to get rid of all the NaN's in the height_2D array produced by the interpolating.
ind = all(~isnan(rastads.height_2D),axis=1)
rasta = rastads.isel(time=np.squeeze(np.where(ind)))

rawp, rastap = xr.align(raw,rasta,join='inner')

tmp=rastap['w_ret']
# tmp=tmp.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
tmp=tmp.where((rastap.Mask_Vz==1))
tmp=tmp.transpose()
fig=plt.figure(figsize=[16,6])    ########### comment me

x = np.tile(rawp.time,(500,1))
y = rastap.height_2D.transpose().values
z = tmp
plt.contourf(x,y,z,np.linspace(-8,8,13),extend='both')

plt.gcf().autofmt_xdate()
plt.gca().xaxis.set_major_formatter(DateFormatter('%H:%M'))

ax=plt.gca()
hcb=plt.colorbar()
hcb.set_label('Updraft m/s')
plt.ylim(0,20)
plt.ylabel('Altitude (km)')
plt.title('Flight 12 updraft Mask 1')

ax2=ax.twinx()
time=rawp.time
iwc=rawp.XKBZR5s
ax2.plot(time,iwc,'k',linewidth=0.3)
ax2.set_ylim(0,10)
_=ax2.set_ylabel('IWC g/m3')


#%%

## Method TWO to integrate different dataset -- with raste-raw-cmb.h5, largest file size, least recommended
## Create updraft - iwc profile

### plt.contourf won't work if the X or Y contains any NaN!!!!
### Therefore need to get rid of all the NaN's in the height_2D array produced by the interpolating.
ind = all(~isnan(cmbds.height_2D),axis=1)
ind = ind & (cmbds.flightnum==13)
ds2 = cmbds.isel(time=np.squeeze(np.where(ind)))
tmp2=ds2['w_ret']
# tmp=tmp.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
## Only mask the variable to be plot, cannot mask all, otherwise height_2D fails to plot contour
tmp2=tmp2.where((ds2.Mask_Vz==1))
# thresh=14
# tmp=tmp.where((tmp<thresh) & (tmp>-thresh))
tmp2=tmp2.transpose()
fig=plt.figure(figsize=[16,6])    ########### comment me

x = np.tile(ds2.time,(500,1))
y = ds2.height_2D.transpose().values
z = tmp2

plt.contourf(x,y,z,np.linspace(-8,8,13),extend='both')
plt.gcf().autofmt_xdate()
plt.gca().xaxis.set_major_formatter(DateFormatter('%H:%M'))

ax=plt.gca()
hcb=plt.colorbar()
hcb.set_label('Updraft m/s')
plt.ylim(0,20)
# plt.xlabel('Time (h)')
plt.ylabel('Altitude (km)')
plt.title('Flight 13 updraft Mask 1')

ax2=ax.twinx()
indflt=(rawvar['flightnum']==13)
iwc=rawvar['XKBZR5s'][indflt]
time=rawvar.index[indflt]
ax2.plot(time,iwc,'k',linewidth=0.3)
ax2.set_ylim(0,10)
_=ax2.set_ylabel('IWC g/m3')
plt.show()


#%%

## Method THREE to integrate different dataset
## Create updraft - iwc profile
ds = loadRastaflt(16)
tmp=ds['w_ret']
# tmp=tmp.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
tmp=tmp.where((ds.Mask_Vz==1))
# thresh=14
# tmp=tmp.where((tmp<thresh) & (tmp>-thresh))
tmp=tmp.transpose()
fig=plt.figure(figsize=[16,6])    ########### comment me
plt.contourf(np.tile(ds.time,(tmp.values.shape[0],1)),ds.height_2D.transpose(),
             tmp,np.linspace(-8,8,13),extend='both')
plt.gcf().autofmt_xdate()
plt.gca().xaxis.set_major_formatter(DateFormatter('%H:%M'))

ax=plt.gca()
hcb=plt.colorbar()
hcb.set_label('Updraft m/s')
plt.ylim(0,20)
# plt.xlabel('Time (h)')
plt.ylabel('Altitude (km)')
plt.title('Flight 16 updraft Mask 1')

ax2=ax.twinx()
indflt=(rawvar['flightnum']==16)
iwc=rawvar['XKBZR5s'][indflt]
time=rawvar.index[indflt]
ax2.plot(time,iwc,'k',linewidth=0.3)
ax2.set_ylim(0,10)
_=ax2.set_ylabel('IWC g/m3')
plt.show()


#%%

plt.scatter(rawvards.air_temperature,rawvards.altitude)
from matplotlib.colors import LogNorm
output = plt.hist2d(rawvards.air_temperature,rawvards.altitude,
                    range=[[-60,0],[200,1400]],bins=np.array([10,12])*5,norm=LogNorm(),cmin=3)
cb = plt.colorbar()
cb.set_label('Frequency')


#%%

attmk=np.array(['no cloud','ice','rain','ice but likely attenuated','ground','ghost ground','interpolated'])
szi=13
for j in range(6):
    szj=j+1
    ds = loadRastaflt(szi)
    tmp=ds.w_ret.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
    tmp=tmp.where(ds.attenuation_phase_flag==szj)
    # tmp=tmp.where((ds.Mask_Vz==1))
    # thresh=14
    # tmp=tmp.where((tmp<thresh) & (tmp>-thresh))
    tmp=tmp.transpose()
    #fig=plt.figure(figsize=[16,6])    ########### comment me
    plt.contourf(np.tile(ds.timeutc,(tmp.values.shape[0],1)),ds.height_2D.transpose(),
                 tmp,np.linspace(-8,8,13),extend='both')
    ax=plt.gca()
    hcb=plt.colorbar()
    hcb.set_label('Updraft m/s')
    plt.ylim(0,20)
    plt.ylabel('Altitude (km)')
    plt.title('Flight '+str(szi)+' Phs Mk '+str(attmk[szj]))
    ax2=ax.twinx()
    indflt=(rawvar['flightnum']==szi)
    iwc=rawvar['XKBZR5s'][indflt]
    time=rawvar.index[indflt]
    ax2.plot(time,iwc,'k',linewidth=0.3)
    ax2.set_ylim(0,10)
    ax2.set_ylabel('IWC g/m3')
    plt.gcf().autofmt_xdate()
    plt.gca().xaxis.set_major_formatter(DateFormatter('%H:%M'))
    plt.show()


#%%

print('Distance between gates is '+str(np.average(np.diff(ds.range)))+' km')


#%%

ds = loadRastaflt(13)

w1=ds.w_wind.values
w2=np.mean( ds.w_ret.values[:,np.array([245,254])], axis=1)

from scipy import stats
plt.scatter(w1,w2,s=.1)
o2o=np.array(plt.gca().get_ylim())*0.5
plt.plot(o2o,o2o,'b')
plt.grid(b=True)
mask = ~np.isnan(w1) & ~np.isnan(w2)
slope, intercept, r_value, p_value, std_err = stats.linregress(w1[mask],w2[mask] )
plt.plot(o2o,o2o*slope+intercept,'r')
plt.xlabel('SAFIRE updraft V5')
plt.ylabel('RASTA updraft')
print('Rsquare '+str(r_value**2))
print('Slope '+str(slope))
print('Intercept '+str(intercept))
plt.legend(['One-to-One','Linear Reg']);


#%%

for i in [12,15]:
    try:
        szi=i+1
        fig=figure(figsize=[16,6])
        ds = loadRastaflt(szi)
        tmp=ds.w_ret.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        #tmp=tmp.where((ds.attenuation_phase_flag==1) | (ds.attenuation_phase_flag==3))
        # tmp=tmp.where((ds.Mask_Vz==1))
        tmp=tmp.transpose()
        plt.contourf(np.tile(ds.timeutc,(tmp.values.shape[0],1)),ds.height_2D.transpose(),
                     tmp,np.linspace(-8,8,13),extend='both')
        ax=plt.gca()
        hcb=plt.colorbar()
        hcb.set_label('Updraft m/s')
        plt.ylim(0,20)
        plt.ylabel('Altitude (km)')
        plt.title('Flight '+str(szi)+' updraft')
        ax2=ax.twinx()
        indflt=(rawvar['flightnum']==szi)
        iwc=rawvar['XKBZR5s'][indflt]
        time=rawvar.index[indflt]
        ax2.plot(time,iwc,'k',linewidth=0.3)
        ax2.set_ylim(0,10)
        ax2.set_ylabel('IWC g/m3')
        plt.gcf().autofmt_xdate()
        plt.gca().xaxis.set_major_formatter(DateFormatter('%H:%M'))
        plt.show()
    except:
        pass


#%%

## vertical updraft profile through time
totaln=23
for i in range(totaln):
    try:
        szi=i+1
        fig=figure(figsize=[16,6])
        ds = loadRastaflt(szi)
#         tmp=ds.w_ret.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        tmp=ds.w_ret.where(ds.Mask_Vz==1)
        tmp=tmp.values.transpose()
        plt.contourf(np.tile(ds.timeutc,(tmp.shape[0],1)),ds.height_2D.values.transpose(),
                     tmp,np.linspace(-8,8,13),extend='both')
        ax=plt.gca()
        hcb=plt.colorbar()
        hcb.set_label('Updraft m/s')
        plt.ylim(0,20)
        plt.ylabel('Altitude (km)')
        plt.title('Flight '+str(szi)+' updraft Mask 1 & 3')
        ax2=ax.twinx()
        indflt=(rawvar['flightnum']==szi)
        iwc=rawvar['XKBZR5s'][indflt]
        time=rawvar.index[indflt]
        ax2.plot(time,iwc,'k',linewidth=0.3)
        ax2.set_ylim(0,10)
        ax2.set_ylabel('IWC g/m3')
        ax2.plot(ax2.get_xlim(),np.array([1.5,1.5]),'k--')
        plt.gcf().autofmt_xdate()
        plt.gca().xaxis.set_major_formatter(DateFormatter('%H:%M'))
        plt.show()
    except:
        pass


#%%

plt.plot(rawvar['altitude'][rawvar['XKBZR5s']>.5].as_matrix())
ylabel('Altitude (m)')


#%%

# grouped w wind profile by MMD cats
totaln=23
cbindbackward=cbmmd=np.zeros([0,1])
cbalt=cbupdf=np.ndarray([500,0])
for i in range(totaln):
# for i in range(12,17):
        szi=i+1
        
        ### For PSD data
        indflt=(raw['flightnum']==szi) & (iwc_mean>1.5)
        tmp = proc['indpsdforward'][indflt]
        tmp = tmp[~np.isnan(tmp)].astype(int)
        #### Need to convert to int otherwise float as index is not acceptable
        mmd = proc['MMD'][ tmp ]
        #### squeeze the dimension otherwise it will create a 3-D array.
        time = raw['Time'][ np.squeeze(proc['indpsdback'][tmp]) ]
        
        ### For RASTA data
        try:
            ds = loadRastaflt(szi).sel(method='nearest',time=np.squeeze(time/3600.),tolerance=1./3600)
        except:
            continue
            pass
#         updf=ds.w_ret.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        if np.prod(ds.w_ret.shape) == 0:
            continue;
        updf = ds.w_ret.where(ds.Mask_Vz==1).values.transpose()
        alt = ds.height_2D.values.transpose()
        if updf.ndim == 1:
            updf=updf[:,np.newaxis]
            alt=alt[:,np.newaxis]
            time=time[np.newaxis,:]
            
        alt = alt-alt[249,:]
        
        cbindbackward=np.concatenate( (cbindbackward,time),axis=0 )
        cbmmd=np.concatenate( (cbmmd,mmd),axis=0 )
        cbupdf=np.concatenate( (cbupdf,updf),axis=1 )
        cbalt=np.concatenate( (cbalt,alt),axis=1 )
#     except:
#         pass


#%%

from mpl_toolkits.basemap import Basemap
m = Basemap(projection='merc',resolution='c',llcrnrlon=110.1,llcrnrlat=-30.1,
            urcrnrlon=150.1,urcrnrlat=5.1)
x,y=m(130,-15)
m.is_land(x,y)


#%%

# land see masks
totaln=23
cbls=cbindbackward=cbmmd=np.zeros([0,1])
cbls=np.zeros(0,)
cbalt=cbupdf=np.ndarray([500,0])
for i in range(totaln):
# for i in range(12,17):
        szi=i+1
        
        ### For PSD data
        indflt=(raw['flightnum']==szi) & (iwc_mean>1.5)
        tmp = proc['indpsdforward'][indflt]
        tmp = tmp[~np.isnan(tmp)].astype(int)
        #### Need to convert to int otherwise float as index is not acceptable
        mmd = proc['MMD'][ tmp ]
        #### squeeze the dimension otherwise it will create a 3-D array.
        time = raw['Time'][ np.squeeze(proc['indpsdback'][tmp]) ]
        lon = raw['lon'][ np.squeeze(proc['indpsdback'][tmp]) ]
        x,y=m(np.squeeze(raw['lon'][ np.squeeze(proc['indpsdback'][tmp]) ]),
              np.squeeze(raw['lat'][ np.squeeze(proc['indpsdback'][tmp]) ]))
        ls=x
        try:
            ttt=len(x)
            for j in range(ttt):
                ls[j]=m.is_land(x[j],y[j])
        except:
            ttt=1
            ls=np.zeros(1,)
            ls[0]=np.array(m.is_land(x,y))
            pass
        
        ### For RASTA data
        try:
            ds = loadRastaflt(szi).sel(method='nearest',time=np.squeeze(time/3600.),tolerance=1./3600)
        except:
            continue
            pass
#         updf=ds.w_ret.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        if np.prod(ds.w_ret.shape) == 0:
            continue;
        updf = ds.w_ret.where(ds.Mask_Vz==1).values.transpose()
        alt = ds.height_2D.values.transpose()
        if updf.ndim == 1:
            updf=updf[:,np.newaxis]
            alt=alt[:,np.newaxis]
            time=time[np.newaxis,:]
#             ls=ls[np.newaxis,:]
            
        alt = alt-alt[249,:]
        
        cbindbackward=np.concatenate( (cbindbackward,time),axis=0 )
        cbmmd=np.concatenate( (cbmmd,mmd),axis=0 )
        cbupdf=np.concatenate( (cbupdf,updf),axis=1 )
        cbalt=np.concatenate( (cbalt,alt),axis=1 )
        cbls=np.concatenate( (cbls,ls),axis=0 )
#     except:
#         pass


#%%

mmdland=cbmmd[cbls==1]
mmdsea=cbmmd[cbls==0]
plt.hist(mmdland)
plt.xlabel('MMD (um)')
plt.ylabel('Frequency')
plt.title('Land PSD')
plt.xlim(0,1200)
plt.show()
plt.hist(mmdsea)
plt.xlabel('MMD (um)')
plt.ylabel('Frequency')
plt.xlim(0,1200)
plt.title('Sea PSD')
plt.show()


#%%

rastap.height_2D[:,249]
# tmpind.shape


#%%

## Method ONE to integrate different dataset -- mostly recommended
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

for flt in [12,13,16]:
    ind = ~isnan(rawvards.mmd)
    ind = ind & (rawvards.flightnum==flt) & (rawvards.XKBZR5s>2.5)# & (rawvards.mmd>650)
    raw = rawvards.isel(time=np.squeeze(np.where(ind)))

    ### plt.contourf won't work if the X or Y contains any NaN!!!!
    ### Therefore need to get rid of all the NaN's in the height_2D array produced by the interpolating.
    ind = all(~isnan(rastads.height_2D),axis=1)
    rasta = rastads.isel(time=np.squeeze(np.where(ind)))

    rawp, rastap = xr.align(raw,rasta,join='inner')

    tmp=rastap['w_ret']
    # tmp=tmp.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
    tmp=tmp.where((rastap.Mask_Vz==1))
    tmp=tmp.transpose()
    fig=plt.figure(figsize=[16,6])    ########### comment me

    mmd = rawp.mmd
    updf = tmp
    alt = rastap.height_2D.where((rastap.Mask_Vz==1)).transpose()
    if flt==16:
        reso = 50
    else:
        reso = 25
    mmdbins = np.percentile(mmd, np.arange(0,100.1,reso))
    mmdind = np.digitize(mmd,mmdbins)

    tmpc=np.array(['r','g','b','k'])
    tmpmmd=np.zeros(len(mmdbins)-1)
    for j in range(1,len(mmdbins)):
        tmpind=np.in1d(mmdind,j)
        x1=updf[:,tmpind]
        y1=alt[:,tmpind]
        x0=np.nanmean(x1,axis=1)
        print(x1.shape)
        y0=np.nanmean(y1,axis=1)
        tmpmmd[j-1]=np.nanmean(mmd[tmpind])
    #             plot(x0,y0,c=colormapping.to_rgba(tmpmmd[j-1]))
        plot(x0,y0,c=tmpc[j-1])

    ax=plt.gca()
    ax.legend(tmpmmd.astype(int).astype(str),loc='best')
    plt.ylabel('Altitude (km)')
    plt.xlabel('Updraft (m/s)')
    plt.xlim(-3,6.5)
    plt.title('Flight '+str(flt)+' updraft IWC cases '+str(len(mmd)))
    fltlvl = nanmean(rastap.height_2D[:,249])
    plt.plot(plt.xlim(),np.array([fltlvl,fltlvl]),'k--')
    plt.plot(np.array([0,0]),plt.ylim(),'k--')
    plt.show()


#%%

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

## PREVIOUS CELL REQUIRED
fig=figure(figsize=[16,6])
mmdbins = np.array([0,200,300,400,500,600,700,800,1000])
mmdind = np.digitize(cbmmd,mmdbins)

colormapping=mpl.cm.ScalarMappable(norm=None,cmap='jet')
colors=colormapping.to_rgba(mmdbins[1:-1])

# tmpc=np.array(['r','g','b','k'])
tmpmmd=np.zeros(len(mmdbins))
for j in range(1,len(mmdbins)-1):## skip the first bin (mmd<200)
    tmpind=np.in1d(mmdind,j+1)
    x1=cbupdf[:,tmpind]
    y1=cbalt[:,tmpind]
    x0=np.nanmean(x1,axis=1)
    y0=np.nanmean(y1,axis=1)
    tmpmmd[j]=np.nanmean(cbmmd[tmpind])
    plot(x0,y0,c=colormapping.to_rgba(tmpmmd[j-1]))
#     plot(x0,y0)

ax=plt.gca()
ax.legend(tmpmmd[1:-1].astype(int).astype(str),loc='best')
plt.ylabel('Elevation (km)')
plt.xlabel('Updraft (m/s)')
plt.xlim(-5,5)
plt.title('All Flights updraft MMD cases '+str(len(np.squeeze(cbmmd))))
plt.plot(plt.xlim(),np.array([0,0]),'k--')
plt.plot(np.array([0,0]),plt.ylim(),'k--')
plt.show()


#%%

### land see masks
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

## PREVIOUS CELL REQUIRED
fig=figure(figsize=[16,6])
mmdbins = np.array([0,200,300,400,500,600,700,800,1000])
mmdind = np.digitize(cbmmd,mmdbins)

colormapping=mpl.cm.ScalarMappable(norm=None,cmap='jet')
colors=colormapping.to_rgba(mmdbins[1:-1])

# tmpc=np.array(['r','g','b','k'])
tmpmmd=np.zeros(len(mmdbins))
for j in range(1,len(mmdbins)-1):## skip the first bin (mmd<200)
    tmpind=np.in1d(mmdind,j+1)
    x1=cbupdf[:,tmpind]
    y1=cbalt[:,tmpind]
    x0=np.nanmean(x1,axis=1)
    y0=np.nanmean(y1,axis=1)
    tmpmmd[j]=np.nanmean(cbmmd[tmpind])
    plot(x0,y0,c=colormapping.to_rgba(tmpmmd[j-1]))
#     plot(x0,y0)

ax=plt.gca()
ax.legend(tmpmmd[1:-1].astype(int).astype(str),loc='best')
plt.ylabel('Elevation (km)')
plt.xlabel('Updraft (m/s)')
plt.xlim(-5,5)
plt.title('All Flights updraft MMD cases '+str(len(np.squeeze(cbmmd))))
plt.plot(plt.xlim(),np.array([0,0]),'k--')
plt.plot(np.array([0,0]),plt.ylim(),'k--')
plt.show()


#%%

## PREVIOUS CELL REQUIRED
fig=figure(figsize=[16,6])
mmdbins = np.array([0,200,300,400,500,600,700,800,1000])
mmdind = np.digitize(cbmmd,mmdbins)

colormapping=mpl.cm.ScalarMappable(norm=None,cmap='jet')
colors=colormapping.to_rgba(mmdbins[1:-1])

# tmpc=np.array(['r','g','b','k'])
tmpmmd=np.zeros(len(mmdbins))
for j in range(1,len(mmdbins)-1):## skip the first bin (mmd<200)
    tmpind=np.in1d(mmdind,j+1)
    x1=cbupdf[:,tmpind]
    y1=cbalt[:,tmpind]
    x0=np.nanstd(x1,axis=1)
    y0=np.nanmean(y1,axis=1)
    tmpmmd[j]=np.nanmean(cbmmd[tmpind])
    plot(x0,y0,c=colormapping.to_rgba(tmpmmd[j-1]))
#     plot(x0,y0)

ax=plt.gca()
ax.legend(tmpmmd[1:-1].astype(int).astype(str),loc='best')
plt.ylabel('Altitude (km)')
plt.xlabel('Updraft (m/s)')
plt.title('All Flights updraft MMD cases '+str(len(np.squeeze(cbmmd))))
plt.plot(plt.xlim(),np.array([0,0]),'k--')
plt.plot(np.array([0,0]),plt.ylim(),'k--')
plt.show()


#%%

## PREVIOUS CELL REQUIRED
from scipy import signal
def moving_average(a, n=3) :
    hlf = int(n/2.)
    n = hlf*2+1
    ret = np.cumsum(a, dtype=float)
    ret[hlf:-hlf]=(ret[n-1:]-ret[:-n+1])/(n-1)
    return np.concatenate((a[:hlf],ret[hlf:-hlf],a[-hlf:]))

def fill_nan(x0):
    ########### Interesting code to interpolate the NaN's
    nans, x= np.isnan(x0), lambda z: z.nonzero()[0]
    x0[nans]= np.interp(x(nans), x(~nans), x0[~nans])
    return x0

fig=figure(figsize=[16,6])
mmdbins = np.array([0,200,300,400,500,600,700,800,1000])
mmdind = np.digitize(cbmmd,mmdbins)

colormapping=mpl.cm.ScalarMappable(norm=None,cmap='jet')
colors=colormapping.to_rgba(mmdbins[1:-1])

# tmpc=np.array(['r','g','b','k'])
tmpmmd=np.zeros(len(mmdbins))
for j in range(1,len(mmdbins)-1):## skip the first bin (mmd<200)
    tmpind=np.in1d(mmdind,j+1)
    x1=cbupdf[:,tmpind]
    y1=cbalt[:,tmpind]
    y0=np.nanmean(y1,axis=1)
    x0=np.nanmean(x1,axis=1)
    
    x0=fill_nan(x0)
    x0=moving_average(x0,n=20)
#     x0=np.diff(x0)
    tmpmmd[j]=np.nanmean(cbmmd[tmpind])
    plot(x0,y0,c=colormapping.to_rgba(tmpmmd[j-1]))
#     plot(x0,y0)

ax=plt.gca()
ax.legend(tmpmmd[1:-1].astype(int).astype(str),loc='best')
plt.ylabel('Altitude (km)')
plt.xlabel('Updraft (m/s)')
# plt.xlim(-.5,.5)
plt.title('All Flights updraft MMD cases '+str(len(np.squeeze(cbmmd))))
plt.plot(plt.xlim(),np.array([0,0]),'k--')
plt.plot(np.array([0,0]),plt.ylim(),'k--')
plt.show()


#%%

# grouped w wind profile by MMD cats
totaln=23
cbindbackward=cbmmd=np.zeros([0,1])
cbalt=cbupdf=np.ndarray([500,0])
# for i in range(totaln):
for i in [11,12,15]:
        szi=i+1
        
        ### For PSD data
        indflt=(raw['flightnum']==szi) & (iwc_mean>1.5)
        tmp = proc['indpsdforward'][indflt]
        tmp = tmp[~np.isnan(tmp)].astype(int)
        #### Need to convert to int otherwise float as index is not acceptable
        mmd = proc['MMD'][ tmp ]
        #### squeeze the dimension otherwise it will create a 3-D array.
        time = raw['Time'][ np.squeeze(proc['indpsdback'][tmp]) ]
        
        ### For RASTA data
        try:
            ds = loadRastaflt(szi).sel(method='nearest',time=np.squeeze(time/3600.),tolerance=1./3600)
        except:
            continue
            pass
#         updf=ds.w_ret.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        if np.prod(ds.w_ret.shape) == 0:
            continue;
        updf = ds.w_ret.where(ds.Mask_Vz==1).values.transpose()
        alt = ds.height_2D.values.transpose()
        if updf.ndim == 1:
            updf=updf[:,np.newaxis]
            alt=alt[:,np.newaxis]
            time=time[np.newaxis,:]
            
        alt = alt-alt[249,:]
        
        fig=figure(figsize=[16,6])
        mmdbins = np.percentile(mmd, np.arange(0,100.1,25))
        mmdind = np.digitize(mmd,mmdbins)

        tmpc=np.array(['r','g','b','k'])
        tmpmmd=np.zeros(len(mmdbins)-1)
        for j in range(1,len(mmdbins)):
            tmpind=np.in1d(mmdind,j)
            x1=updf[:,tmpind]
            y1=alt[:,tmpind]
            x0=np.nanmean(x1,axis=1)
            y0=np.nanmean(y1,axis=1)
            tmpmmd[j-1]=np.nanmean(mmd[tmpind])
        #             plot(x0,y0,c=colormapping.to_rgba(tmpmmd[j-1]))
            plot(x0,y0,c=tmpc[j-1])

        ax=plt.gca()
        ax.legend(tmpmmd.astype(int).astype(str),loc='best')
        plt.ylabel('Altitude (km)')
        plt.xlabel('Updraft (m/s)')
        plt.xlim(-5,5)
        plt.title('Flight '+str(szi)+' updraft IWC cases '+str(len(mmd)))
        plt.plot(plt.xlim(),np.array([0,0]),'k--')
        plt.plot(np.array([0,0]),plt.ylim(),'k--')
        plt.show()


#%%

# grouped w wind profile by MMD cats
totaln=23
cbindbackward=cbmmd=np.zeros([0,1])
cbalt=cbupdf=np.ndarray([500,0])
# for i in range(totaln):
for i in [11,12,15]:
        szi=i+1
        
        ### For PSD data
        indflt=(raw['flightnum']==szi) & (iwc_mean>1.5)
        tmp = proc['indpsdforward'][indflt]
        tmp = tmp[~np.isnan(tmp)].astype(int)
        #### Need to convert to int otherwise float as index is not acceptable
        mmd = proc['MMD'][ tmp ]
        #### squeeze the dimension otherwise it will create a 3-D array.
        time = raw['Time'][ np.squeeze(proc['indpsdback'][tmp]) ]
        
        ### For RASTA data
        try:
            ds = loadRastaflt(szi).sel(method='nearest',time=np.squeeze(time/3600.),tolerance=1./3600)
        except:
            continue
            pass
#         updf=ds.w_ret.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        if np.prod(ds.w_ret.shape) == 0:
            continue;
        updf = ds.w_ret.where(ds.Mask_Vz==1).values.transpose()
        alt = ds.height_2D.values.transpose()
        if updf.ndim == 1:
            updf=updf[:,np.newaxis]
            alt=alt[:,np.newaxis]
            time=time[np.newaxis,:]
            
#         alt = alt-alt[249,:]
        
        fig=figure(figsize=[16,6])
        mmdbins = np.percentile(mmd, np.arange(0,100.1,25))
        mmdind = np.digitize(mmd,mmdbins)

        tmpc=np.array(['r','g','b','k'])
        tmpmmd=np.zeros(len(mmdbins)-1)
        for j in range(1,len(mmdbins)):
            tmpind=np.in1d(mmdind,j)
            x1=updf[:,tmpind]
            y1=alt[:,tmpind]
            x0=np.nanmean(x1,axis=1)
            y0=np.nanmean(y1,axis=1)
            tmpmmd[j-1]=np.nanmean(mmd[tmpind])
        #             plot(x0,y0,c=colormapping.to_rgba(tmpmmd[j-1]))
            plot(x0,y0,c=tmpc[j-1])

        ax=plt.gca()
        ax.legend(tmpmmd.astype(int).astype(str),loc='best')
        plt.ylabel('Altitude (km)')
        plt.xlabel('Updraft (m/s)')
        plt.xlim(-5,5)
        plt.title('Flight '+str(szi)+' updraft IWC cases '+str(len(mmd)))
#         plt.plot(plt.xlim(),np.array([0,0])+12,'k--')
        plt.plot(np.array([0,0]),plt.ylim(),'k--')
        plt.show()


#%%

# grouped w wind profile by IWC cats
import pandas as pd


totaln=23
for i in range(totaln):
# for i in range(12,13):
    try:
        szi=i+1
        fig=figure(figsize=[16,6])
        
        ### For PSD data
        indflt=((raw['flightnum']==szi) & (iwc_mean>1.5))
        #### Need to convert to int otherwise float as index is not acceptable
        mmd = iwc_mean[ indflt ]
        #### squeeze the dimension otherwise it will create a 3-D array.
        time = raw['Time'][ indflt ]
        
        ### For RASTA data
        ds = loadRastaflt(szi).sel(method='nearest',time=np.squeeze(time/3600.),tolerance=1./3600)
#         updf=ds.w_ret.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        updf = ds.w_ret.where(ds.Mask_Vz==1)
        
        mmdbins = np.percentile(mmd, np.arange(0,100.1,25))
        mmdind = np.digitize(mmd,mmdbins)
        colormapping=mpl.cm.ScalarMappable(norm=None, cmap='jet')
        colors=colormapping.to_rgba(mmdbins[:-1]+np.diff(mmdbins))
        
        tmpc=np.array(['r','g','b','k'])
        tmpmmd=np.zeros(len(mmdbins)-1)
        for j in range(1,len(mmdbins)):
            tmpind=np.in1d(mmdind,j)
            x1=updf.values[tmpind,:].transpose()
            tmp=ds.height_2D.values
            y1=(tmp[tmpind,:]-tmp[tmpind,249][:,np.newaxis]).transpose()
            x0=np.nanmean(x1,axis=1)
            y0=np.nanmean(y1,axis=1)
            tmpmmd[j-1]=np.nanmean(mmd[tmpind])
#             plot(x0,y0,c=colormapping.to_rgba(tmpmmd[j-1]))
            plot(x0,y0,c=tmpc[j-1])
            
        ax=plt.gca()
        ax.legend(tmpmmd.astype(str),loc='best')
        plt.ylabel('Altitude (km)')
        plt.xlabel('Updraft (m/s)')
        plt.title('Flight '+str(szi)+' updraft IWC cases '+str(len(np.squeeze(mmd))))
        plt.plot(plt.xlim(),np.array([0,0]),'k--')
        plt.plot(np.array([0,0]),plt.ylim(),'k--')
        plt.show()
    except:
        pass

