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

get_ipython().run_cell_magic('capture', '', "%pylab inline\nimport xarray as xr\nimport glob\n\ndef loadRastaflt(szi):\n#     matlabpath='/data/mcfarq/a/szhu28/research/HIWC/10_150901Wholeflight/src/analysis/example/data/'\n    rastafn=glob.glob(datapath+'*_F'+str(szi)+'_*.nc')\n    return xr.open_dataset(rastafn[0])")


#%%

import scipy.io as sio    # For .mat version before 7.3
import h5py    # For .mat version after 7.3
import pandas as pd
matlabpath='/data/mcfarq/a/szhu28/research/HIWC/10_150901Wholeflight/src/analysis/example/data/'
datapath='/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/RASTA/data/'
filenames=np.array(['processed.mat','rawFIT.mat','rawPSD.mat','rawVAR.mat'])
'''
for szi in range(4):
    tmpfn=matlabpath+filenames[szi]
    try:
        raw=h5py.File(tmpfn,'r')
#         print('hdf '+str(szi))
        raw.close()
    except:
        raw=sio.loadmat(tmpfn)
#         raw=sio.whosmat(tmpfn)
#         varnames=raw.items()
'''
raw=sio.loadmat(matlabpath+filenames[3])
proc=sio.loadmat(matlabpath+filenames[0])
proc['indpsdforward']-=1                     ############ PYTHON VS MATLAB
proc['indpsdback']-=1

### 5 seconds running mean of iwc
tmpiwc=pd.DataFrame(raw['TWCIKPZRgm3'].ravel())
iwc_mean=tmpiwc.rolling(window=5,center=True,min_periods=3).mean().as_matrix()


#%%

ds = loadRastaflt(13)
tmp=ds['Vz']
# tmp=tmp.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
tmp=tmp.where((ds.Mask_Vz==1))
# thresh=14
# tmp=tmp.where((tmp<thresh) & (tmp>-thresh))
tmp=tmp.transpose()
fig=plt.figure(figsize=[16,6])    ########### comment me
plt.contourf(np.tile(ds.time.values,(tmp.values.shape[0],1)),ds.height_2D.transpose(),
             tmp,np.linspace(-8,8,13),extend='both')
ax=plt.gca()
hcb=plt.colorbar()
hcb.set_label('Updraft m/s')
plt.ylim(0,20)
plt.xlabel('Time (h)')
plt.ylabel('Altitude (km)')
plt.title('Flight 13 updraft Mask 1')
ax2=ax.twinx()
indflt=(raw['flightnum']==13)
iwc=raw['TWCIKPZRgm3'][indflt]
time=raw['Time'][indflt]
ax2.plot(time/3600.,iwc,'k',linewidth=0.3)
ax2.set_ylim(0,10)
_=ax2.set_ylabel('IWC g/m3')


#%%

attmk=np.array(['no cloud','ice','rain','ice but likely attenuated','ground','ghost ground','interpolated'])
szi=13
for j in range(6):
    szj=j+1
    ds = loadRastaflt(szi)
    tmp=ds.Vz.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
    tmp=tmp.where(ds.attenuation_phase_flag==szj)
    # tmp=tmp.where((ds.Mask_Vz==1))
    # thresh=14
    # tmp=tmp.where((tmp<thresh) & (tmp>-thresh))
    tmp=tmp.transpose()
    #fig=plt.figure(figsize=[16,6])    ########### comment me
    plt.contourf(np.tile(ds.time.values,(tmp.values.shape[0],1)),ds.height_2D.transpose(),
                 tmp,np.linspace(-8,8,13),extend='both')
    ax=plt.gca()
    hcb=plt.colorbar()
    hcb.set_label('Updraft m/s')
    plt.ylim(0,20)
    plt.xlabel('Time (h)')
    plt.ylabel('Altitude (km)')
    plt.title('Flight '+str(szi)+' Phs Mk '+str(attmk[szj]))
    ax2=ax.twinx()
    indflt=(raw['flightnum']==szi)
    iwc=raw['TWCIKPZRgm3'][indflt]
    time=raw['Time'][indflt]
    ax2.plot(time/3600.,iwc,'k',linewidth=0.3)
    ax2.set_ylim(0,10)
    ax2.set_ylabel('IWC g/m3')
    plt.show()


#%%

tmp=(~np.isnan(ds['Vz'].values)).sum(axis=0)
plt.step(range(len(tmp)),tmp)
plt.xlabel('Vertical Gates')
plt.ylabel('Effective data points')
print('Gate    '+str(range(245,255)))
print('Points '+str(tmp[245:255]))


#%%

for i in range(1):
    szi=i+1
    try:
        ds=loadRastaflt(szi)
    except:
        continue;
    tmp=ds['height_2D']-np.mean( ds['height_2D'][:,np.array([249,250])],axis=1)
    tmp=tmp.transpose()
    plt.contour(np.tile(ds.time.values,(tmp.values.shape[0],1)),
                np.tile(np.linspace(0,499,500)[:,np.newaxis],(1,tmp.values.shape[1])) ,tmp,
                np.linspace(-12,12,7),extend='both')
    plt.ylabel('Vertical Gates')
    plt.xlabel('Time (seconds)')
    plt.title('Flt '+str(szi)+' Gates relative height')

    ax=plt.gca()
    hcb=plt.colorbar()
    hcb.set_label('Gate relative Altitude (m)')
    plt.xlabel('Time (h)')
    plt.ylabel('Gate number')
    ax2=ax.twinx()
    indflt=(raw['flightnum']==szi)
    iwc=raw['roll'][indflt]
    time=raw['Time'][indflt]
    ax2.plot(time/3600.,iwc,'k',linewidth=0.3)
    #ax2.set_ylim(0,10)
    ax2.set_ylabel('Aircraft Roll angle')
    plt.show()


#%%

print('Distance between gates is '+str(np.average(np.diff(ds.range)))+' km')


#%%

for j in range(23):
    try:
        szj=j+1
        ds = loadRastaflt(szj)
        indflt=(raw['flightnum']==szj)
        wwind=raw['wwind'][indflt]
        time=raw['Time'][indflt]

        timera=ds.time.values*3600
        timerauni=np.unique(timera.round())
        timerauni.shape

        common=numpy.intersect1d(time, timerauni)
        ind1=np.array([ np.sum(a == common) for a in time ]).astype(bool)
        ind2=np.array([ np.sum(a == common) for a in timerauni ]).astype(bool)
        w1=wwind[ind1]
        w2=ds.w_wind.values[ind2]

        from scipy import stats
        plt.scatter(w1,w2,s=.1)
        o2o=np.array(plt.gca().get_ylim())*0.5
        plt.plot(o2o,o2o,'b')
        plt.grid(b=True)
        mask = ~np.isnan(w1) & ~np.isnan(w2)
        slope, intercept, r_value, p_value, std_err = stats.linregress(w1[mask],w2[mask] )
        plt.plot(o2o,o2o*slope+intercept,'r')
        plt.xlabel('SAFIRE updraft V3')
        plt.ylabel('SAFIRE updraft V5')
        print('Flight '+str(szj))
        print('\tRsquare '+str(r_value**2))
        print('\t\tSlope '+str(slope))
        print('\t\t\tIntercept '+str(intercept))
        plt.legend(['One-to-One','Linear Reg']);
        if szj!=23:
            plt.close()
    except:
        pass


#%%

ds = loadRastaflt(13)
indflt=(raw['flightnum']==13)
wwind=raw['wwind'][indflt]
time=raw['Time'][indflt]

timera=ds.time.values*3600
timerauni=np.unique(timera.round())
timerauni.shape

common=numpy.intersect1d(time, timerauni)
ind1=np.array([ np.sum(a == common) for a in time ]).astype(bool)
ind2=np.array([ np.sum(a == common) for a in timerauni ]).astype(bool)
w1=wwind[ind1]
w2=np.mean( ds.Vz.values[ind2][:,np.array([245,254])], axis=1)

from scipy import stats
plt.scatter(w1,w2,s=.1)
o2o=np.array(plt.gca().get_ylim())*0.5
plt.plot(o2o,o2o,'b')
plt.grid(b=True)
mask = ~np.isnan(w1) & ~np.isnan(w2)
slope, intercept, r_value, p_value, std_err = stats.linregress(w1[mask],w2[mask] )
plt.plot(o2o,o2o*slope+intercept,'r')
plt.xlabel('SAFIRE updraft V3')
plt.ylabel('RASTA updraft')
print('Rsquare '+str(r_value**2))
print('Slope '+str(slope))
print('Intercept '+str(intercept))
plt.legend(['One-to-One','Linear Reg']);


#%%

ds = loadRastaflt(13)

w1=ds.w_wind.values
w2=np.mean( ds.Vz.values[:,np.array([245,254])], axis=1)

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
        tmp=ds.Vz.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        #tmp=tmp.where((ds.attenuation_phase_flag==1) | (ds.attenuation_phase_flag==3))
        # tmp=tmp.where((ds.Mask_Vz==1))
        tmp=tmp.transpose()
        plt.contourf(np.tile(ds.time.values,(tmp.values.shape[0],1)),ds.height_2D.transpose(),
                     tmp,np.linspace(-8,8,13),extend='both')
        ax=plt.gca()
        hcb=plt.colorbar()
        hcb.set_label('Updraft m/s')
        plt.ylim(0,20)
        plt.xlabel('Time (h)')
        plt.ylabel('Altitude (km)')
        plt.title('Flight '+str(szi)+' updraft')
        ax2=ax.twinx()
        indflt=(raw['flightnum']==szi)
        iwc=raw['TWCIKPZRgm3'][indflt]
        time=raw['Time'][indflt]
        ax2.plot(time/3600.,iwc,'k',linewidth=0.3)
        ax2.set_ylim(0,10)
        ax2.set_ylabel('IWC g/m3')
        plt.show()
    except:
        pass


#%%

totaln=23
for i in range(totaln):
    try:
        szi=i+1
        fig=figure(figsize=[16,6])
        ds = loadRastaflt(szi)
#         tmp=ds.Vz.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        tmp=ds.Vz.where(ds.Mask_Vz==1)
        tmp=tmp.values.transpose()
        plt.contourf(np.tile(ds.time.values,(tmp.shape[0],1)),ds.height_2D.values.transpose(),
                     tmp,np.linspace(-8,8,13),extend='both')
        ax=plt.gca()
        hcb=plt.colorbar()
        hcb.set_label('Updraft m/s')
        plt.ylim(0,20)
        plt.xlabel('Time (h)')
        plt.ylabel('Altitude (km)')
        plt.title('Flight '+str(szi)+' updraft Mask 1 & 3')
        ax2=ax.twinx()
        indflt=(raw['flightnum']==szi)
        iwc=raw['TWCIKPZRgm3'][indflt]
        time=raw['Time'][indflt]
        ax2.plot(time/3600.,iwc,'k',linewidth=0.3)
        ax2.set_ylim(0,10)
        ax2.set_ylabel('IWC g/m3')
        ax2.plot(ax2.get_xlim(),np.array([1.5,1.5]),'k--')
        plt.show()
    except:
        pass


#%%

# w wind profile tripcolor plot
totaln=23
cbx=cby=cbz=np.array([])
for i in range(totaln):
# for i in range(12,13):
    try:
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
        ds = loadRastaflt(szi).sel(method='nearest',time=np.squeeze(time/3600.),tolerance=1./3600)
#         tmp=ds.Vz.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        updf=ds.Vz.where(ds.Mask_Vz==1)
        
#         ### For some reason we don't need to provide norm=cnorm for it to normalize mmd
#         ### But once the colorv.to_rgba(mmd) is initilized by mmd, it no longer chanegs and doesn't
#         ### give the correct scale for mmd/2.
#         # Link http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.scatter
#         # cnorm=mpl.colors.Normalize()
#         # cnorm.autoscale(mmd.ravel())
#         colorv=mpl.cm.ScalarMappable(norm=None, cmap=cmap)
#         colorv.to_rgba(mmd);
#         # colorv.to_rgba(mmd/2)
        
        x1=updf.transpose().values
        y1=ds.height_2D.transpose().values
        y1=y1-y1[249,:]
        z1=np.tile(squeeze(mmd),[x1.shape[0],1])
        x=x1.ravel()
        y=y1.ravel()
        z=z1.ravel()
        indvalid=(~np.isnan(x)) & ~np.isnan(y) & ~np.isnan(z)
        cbx=np.concatenate((cbx,x[indvalid]))
        cby=np.concatenate((cby,y[indvalid]))
        cbz=np.concatenate((cbz,z[indvalid]))
    except:
        pass


#%%

plot(raw['altGPS'][iwc_mean>1.5])
ylabel('Altitude (m)')


#%%

fig=figure(figsize=[16,6])
cmap=mpl.cm.get_cmap()
plt.tricontourf(cbx,cby,cbz,64,cmap=cmap,vmin=250,vmax=850,extend='both')
hcb=plt.colorbar()
hcb.set_label('MMD (um)')
ax=plt.gca()
plt.ylabel('Altitude (km)')
plt.xlabel('Updraft (m/s)')
plt.title('All Flights updraft profile cat by MMD')
plt.plot(plt.xlim(),np.array([0,0]),'k--')
plt.plot(np.array([0,0]),plt.ylim(),'k--')
plt.grid(b=True)
plt.show()


#%%

# grouped w wind profile by MMD cats
totaln=23
cbindbackward=cbmmd=np.zeros([0,1])
cbalt=cbupdf=np.ndarray([500,0])
for i in range(totaln):
# for i in range(12,13):
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
#         updf=ds.Vz.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        if np.prod(ds.Vz.shape) == 0:
            continue;
        updf = ds.Vz.where(ds.Mask_Vz==1).values.transpose()
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

fig=figure(figsize=[16,6])
mmdbins = np.percentile(cbmmd, np.arange(0,100.1,25))
mmdind = np.digitize(cbmmd,mmdbins)

tmpc=np.array(['r','g','b','k'])
tmpmmd=np.zeros(len(mmdbins)-1)
for j in range(1,len(mmdbins)):
    tmpind=np.in1d(mmdind,j)
    x1=cbupdf[:,tmpind]
    y1=cbalt[:,tmpind]
    x0=np.nanmean(x1,axis=1)
    y0=np.nanmean(y1,axis=1)
    tmpmmd[j-1]=np.nanmean(cbmmd[tmpind])
#             plot(x0,y0,c=colormapping.to_rgba(tmpmmd[j-1]))
    plot(x0,y0,c=tmpc[j-1])

ax=plt.gca()
ax.legend(tmpmmd.astype(int).astype(str),loc='best')
plt.ylabel('Altitude (km)')
plt.xlabel('Updraft (m/s)')
plt.title('All Flights updraft MMD cases '+str(len(np.squeeze(cbmmd))))
plt.plot(plt.xlim(),np.array([0,0]),'k--')
plt.plot(np.array([0,0]),plt.ylim(),'k--')
plt.show()


#%%

sum(indflt)


#%%

### Vertical wind profiles in different IWC categories

# grouped w wind profile by MMD cats
import numpy.core.defchararray as npch

mmdbins = np.array([0,300,400,500,600,700,800,1000,2000])
mmdind = np.digitize(mmd,mmdbins)
colormapping=mpl.cm.ScalarMappable(norm=None, cmap='jet')
colors=colormapping.to_rgba( np.arange(len(mmdbins)-1) )

totaln=23
for i in range(totaln):
# for i in range(12,13):
#     try:
        szi=i+1
        fig=figure(figsize=[16,6])
        
        ### For PSD data
        indflt=(raw['flightnum']==szi) & (iwc_mean>1.5)
        tmp = proc['indpsdforward'][indflt]
        tmp = tmp[~np.isnan(tmp)].astype(int)
        #### Need to convert to int otherwise float as index is not acceptable
        mmd = proc['MMD'][ tmp ]
        #### squeeze the dimension otherwise it will create a 3-D array.
        time = raw['Time'][ np.squeeze(proc['indpsdback'][tmp]) ]
        
        ### For RASTA data
        ds = loadRastaflt(szi).sel(method='nearest',time=np.squeeze(time/3600.),tolerance=1./3600)
        tmp=ds.height_2D.values
#         updf=ds.Vz.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        updf = ds.Vz.where(ds.Mask_Vz==1)
        
#         mmdbins = np.percentile(mmd, np.arange(0,100.1,25))
        
        for j in range(1,len(mmdbins)):
            tmpind=np.in1d(mmdind,j)
            if np.sum(tmpind) ==0:
                plot([],[],c=colors[j-1])
                continue
            x1=updf.values[tmpind,:].transpose()
            y1=(tmp[tmpind,:]-tmp[tmpind,249][:,np.newaxis]).transpose()
            x0=np.nanmean(x1,axis=1)
            y0=np.nanmean(y1,axis=1)
#             plot(x0,y0,c=colormapping.to_rgba(tmpmmd[j-1]))
            plot(x0,y0,c=colors[j-1])
            
        ax=plt.gca()
        ax.legend( npch.add(npch.add( mmdbins[:-1].astype(str),' - '),mmdbins[1:].astype(str) ),loc='best')
        plt.ylabel('Altitude (km)')
        plt.xlabel('Updraft (m/s)')
        plt.title('Flight '+str(szi)+' updraft MMD cases '+str(len(np.squeeze(mmd))))
        plt.plot(plt.xlim(),np.array([0,0]),'k--')
        plt.plot(np.array([0,0]),plt.ylim(),'k--')
        plt.show()
#     except:
#         pass


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
#         updf=ds.Vz.where((ds.Mask_Vz==1) | (ds.Mask_Vz==3))
        updf = ds.Vz.where(ds.Mask_Vz==1)
        
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


#%%

mmd.shape


#%%

totaln=23
for i in range(totaln):
    szi=i+1
    indflt = (iwc_mean>1.5) & (raw['flightnum']==i+1) & (~np.isnan(proc['indpsdforward']))
    time=raw['Time'][indflt]
    try:
        ds = loadRastaflt(szi).w_wind.sel(method='nearest',time=np.squeeze(time/3600.),tolerance=1./3600)
    except:
        pass
    tmp=ds.values
    try:
        if tmp.shape[0]==time.shape[0]:
            pass
        else:
            print('error in flt '+str(szi))
            continue
    except:
        pass
        continue
    
    wwind = tmp
    
    tmp = proc['indpsdforward'][indflt]
    tmp = tmp[~np.isnan(tmp)].astype(int)
    #### Need to convert to int otherwise float as index is not acceptable
    mmd = proc['MMD'][ tmp ]
    iwc = iwc_mean[ proc['indpsdback'][tmp] ]
    if np.prod(mmd.shape)==0:
        print('0 mmd valid in flt '+str(szi))
        continue
        
    scatter(mmd,wwind,c=iwc,s=5,lw=0)
    hcb=plt.colorbar()
    hcb.set_label('IWC (g/m3)')
    plt.ylabel('Rasta V5 w wind (m/s)')
    plt.xlabel('MMD (um)')
    ax=plt.gca()
#     ax.set_xlim(1,5)
    ax.set_ylim(-15,15)
    plt.title('Flight '+str(szi)+' HWIC cases '+str(len(np.squeeze(mmd))))
    plt.show()


#%%

tmp

