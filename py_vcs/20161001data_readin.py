##%%writefile jupyter_notebook_exporter_python_mycfg.tpl

# coding: utf-8

#%%

get_ipython().run_line_magic('pylab', 'inline')
import pandas as pd
import xarray as xr
import glob, re, os
import scipy
from scipy import io, optimize
import pickle
import datetime

## For debug mode
from IPython.core.debugger import Tracer
#Tracer()() #this one triggers the debugger


#%%

numofflts = 23
rootpath = '/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin'


#%%

## For IKP2 V5 data
def read_ikp2(filename):

    tmp = pd.read_csv(filename,header=4,usecols=range(14),
                      na_values={'SIAltm':99999.,
                                 'SINShead':999.,
                                 'SRHWVSS':999.,
                                 'Swdir':999.,
                                 'Swspd':999.,
                                 'XKBZR5s':-999.})
    tmptime = pd.to_timedelta(tmp['Stimech'])
    tmptime = pd.to_datetime(datestr) + tmptime
    tmp['Stimech'] = tmptime

    tmp.set_index(['Stimech'],verify_integrity=False,inplace=True,drop=False)
    tmp.index.name=None
    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in IKP is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]
    return tmp

## For SAFFIRE V4 data
def read_saffire(filename):

    '''
    Things to note
    1. The time is in seconds and can be non-sharp seconds.
    '''
    
    ##     Tracer()() #this one triggers the debugger
    varnames = [
        'Timeinsecond',
        'event_marker',
        'latitude',
        'longitude',
        'altitude',
        'altitude',
        'platform_roll_angle',
        'platform_pitch_angle',
        'platform_orientation',
        'air_pressure',
        'air_temperature',
        'air_temperature',
        'air_temperature',
        'dew_point_temperature',
        'relative_humidity',
        'humidity_mixing_ratio',
        'humidity_mixing_ratio',
        'humidity_mixing_ratio',
        'platform_speed_wrt_air',
        'platform_acceleration',
        'platform_course',
        'platform_speed_wrt_ground',
        'platform_course',
        'platform_speed_wrt_ground',
        'upward_platform_speed_wrt_ground',
        'angle_of_attack',
        'angle_of_sideslip',
        'eastward_wind',
        'northward_wind',
        'upward_air_velocity',
        'wind_from_direction',
        'wind_speed',
        'mic_msofreqice_rs_sync_1']

    lookup = 'Warning : most measurements are not valid before take-off and after landing'
    comments = []
    with open(filename) as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup in line:
                skipline = num
                break
            if num>50:
                comments.append(line)
                
    [comments.pop() for i in range(4)]
    comments = ''.join(comments)
    ### r"\s+" refers to one or more occurences of whitespace, while r"\s*" will match zero and would raise a warning
    tmp = pd.read_csv(filename,skiprows=skipline,names=varnames,sep=r"\s+",na_values=3.40282347e+38)
    tmp['timeutc'] = pd.to_datetime(
        pd.Series(np.array(round(tmp['Timeinsecond']), dtype='timedelta64[s]')))+\
        (pd.to_datetime(datestr)-pd.to_datetime('1970-01-01'))

    tmp.set_index(['timeutc'],verify_integrity=False,inplace=True,drop=True)
    tmp.index.name=None
    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in saffire is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]

    return tmp, comments

## For Robust data
def read_robust(filename):
    ## Use explicit names, since the name in data file could  be inconsistent
    varnames = ['GMT','TWC_robust']
    tmp = pd.read_excel(filename,names=varnames)
    tmp['timeutc'] = pd.DatetimeIndex(tmp['GMT']).round('1s')
    tmp.set_index(['timeutc'],inplace=True,drop=True,verify_integrity=False)
    tmp.index.name=None
    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in Robust is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]
    return tmp

## For PSD / MSD data
def read_sd(filename):
    
    tmp = pd.read_csv(filename,header=0,sep=r"\s+")
    tmptime = pd.to_timedelta(tmp['-999'],unit='s')
    tmptime = pd.to_datetime(datestr) + tmptime
    tmp['-999'] = tmptime
    tmp.set_index(['-999'],verify_integrity=False,inplace=True,drop=True)
    tmp.index.name=None
    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in IKP is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]
    return tmp

## For counts data
def read_count(filename):
    
    tmp = pd.read_csv(filename,header=None)
    ### Use tmp[0] instead of tmp['0'] here because index type here is integer rather than string
    tmptime = pd.to_timedelta(tmp[0],unit='s')
    tmptime = pd.to_datetime(datestr) + tmptime
    tmp[0] = tmptime
    tmp.set_index([0],verify_integrity=False,inplace=True,drop=True)
    tmp.index.name=None
    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in IKP is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]
    return tmp


#%%

## Create meta file information
### For flight dates and SAFFIRE V4 files
datapath = rootpath + '/data_src/data20161012/saffire'
flist = np.array( sorted(glob.glob(datapath+'/F20*.txt')) )

datestrs = pd.Series(np.empty(numofflts))
saffirefn = pd.Series(np.empty(numofflts))

for szi in range(numofflts):
    fn = os.path.basename(flist[szi])
    saffirefn[szi] = fn
    tmpdate = re.search(r'(?<=v4_).*(?=_)',fn).group()
    tmpdate = tmpdate[:4]+'-'+tmpdate[4:6]+'-'+tmpdate[6:]
    datestrs[szi] = tmpdate

fileinfo = pd.DataFrame(index=np.arange(1,24).astype(str))
fileinfo['flightdate'] = np.array(datestrs)
fileinfo['comments'] = NaN
fileinfo['saffire'] = np.array(saffirefn)

### For IKP V5 files
datapath = rootpath + '/data_src/data20161012/ikp'
flist = sorted( glob.glob(datapath+'/f20*.csv') )
flist = [os.path.basename(x) for x in flist]
flist.insert(4,None) # no data of ikp from flight 5
fileinfo['ikp'] = flist

### For Robust TWC
datapath = rootpath + '/data_src/data20151203/dataHAIC'
flist =  sorted(glob.glob(datapath+'/Robust_data_flt*.xls'))
flist = [os.path.basename(x) for x in flist]
fileinfo['robust'] = flist

### For PSD

### These 2 variables are used in the next cell as well.
sdnamelist = ['psddmax','psddeq','msddmax','msddeq','psdly']
countnamelist = ['2dscounts','pipcounts']

datapath = rootpath + '/data_src/data20151203/dataHAIC'

tmp = [os.path.basename(x) for x in sorted(glob.glob(datapath+r'/*F#*Composite_Dmax*.txt'))]
tmp.insert(20,None) # no psd data from flight 21
fileinfo[sdnamelist[0]] = tmp

tmp = [os.path.basename(x) for x in sorted(glob.glob(datapath+r'/*F#*Composite_Deq*.txt'))]
tmp.insert(20,None)
fileinfo[sdnamelist[1]] = tmp

tmp = [os.path.basename(x) for x in sorted(glob.glob(datapath+r'/*F#*MassSizeD_Dmax*.txt'))]
tmp.insert(20,None)
fileinfo[sdnamelist[2]] = tmp

tmp = [os.path.basename(x) for x in sorted(glob.glob(datapath+r'/*F#*MassSizeD_Deq*.txt'))]
tmp.insert(20,None)
fileinfo[sdnamelist[3]] = tmp

tmp = [os.path.basename(x) for x in sorted(glob.glob(datapath+r'/*F#*Composite_Ly*.txt'))]
tmp.insert(20,None)
fileinfo[sdnamelist[4]] = tmp

### For counts
datapath = rootpath + '/../counts_HIWC/Darwin/Darwin'
tmp = [os.path.basename(x) for x in sorted(glob.glob(datapath+r'/*chx0_v*.csv'))]
[tmp.insert(szi-1,None) for szi in [14,21,22] ]
fileinfo[countnamelist[0]] = tmp
tmp = [os.path.basename(x) for x in sorted(glob.glob(datapath+r'/*Comptage*.csv'))]
tmp.insert(20,None)
fileinfo[countnamelist[1]] = tmp


#%%

## Create rawvar and rawpsd variables
rawvar = pd.DataFrame()
rawpsd = pd.DataFrame()
for flt in range(1,numofflts+1):
# for flt in range(1,1+1):
    print('Processing flight '+str(flt)+'...')
    filenames = fileinfo.loc[str(flt),:]
    datestr = filenames['flightdate'] # Global function that will be used in the read file functions
    ### IKP
    datapath = rootpath + '/data_src/data20161012/ikp'
    filename = datapath + '/' + str(filenames['ikp'])
    if os.path.isfile(filename):
        tmp1 = read_ikp2(filename)
    else:
        tmp1 = None

    ### SAFFIRE
    datapath = rootpath + '/data_src/data20161012/saffire'
    filename = datapath + '/' + str(filenames['saffire'])
    if os.path.isfile(filename):
        tmp2, comments = read_saffire(filename)
        fileinfo.loc[str(flt),'comments'] = comments
    else:
        tmp2 = None

    ### Robust
    datapath = rootpath + '/data_src/data20151203/dataHAIC'
    filename = datapath + '/' + str(filenames['robust'])
    if os.path.isfile(filename):
        tmp3 = read_robust(filename)
        fileinfo.loc[str(flt),'comments'] = comments
    else:
        tmp3 = None

    ### SD
    datapath = rootpath + '/data_src/data20151203/dataHAIC'
    tmpsd = []
    for varname in sdnamelist:
        filename = datapath + '/' + str(filenames[varname])
        if os.path.isfile(filename):
            tmpsd.append(read_sd( filename ))
        else:
            print('Missing PSD file')
            tmpsd.append(None)
    ### counts
    ### counts are omitted for now because it's in 1-sec resolution. Need future work.
    # datapath = rootpath + '/../counts_HIWC/Darwin/Darwin'
    
    result = pd.concat([tmp1, tmp2, tmp3], axis=1,verify_integrity=True)
    ### Drop the first line where all data is NaT or NaN, this must be done before creating indpsdbackward and flightnum
    result.dropna(axis=0,how='all',inplace=True)
    result['flightnum'] = flt
    rawvar = pd.concat([rawvar,result],verify_integrity=True,axis=0)
    if all([x is None for x in tmpsd]):
        sdresult = None
    else:
        sdresult = pd.concat(tmpsd, axis=1,verify_integrity=True,keys=sdnamelist)
    rawpsd = pd.concat([rawpsd,sdresult], axis=0,verify_integrity=True)
    
### combine sd after all flights are processed!
### Create indpsdforward
tmpind = rawpsd.shape
rawvar = pd.concat([rawvar,pd.Series(np.arange(tmpind[0]),index=rawpsd.index,
                                         name='indpsdforward')],axis=1,verify_integrity=True)

### Create indpsdbackward
tmpshp = rawvar.shape
tmpshp = tmpshp[0]
tmpseries = np.arange(tmpshp)
indpsdforward = rawvar['indpsdforward'].as_matrix()
indpsdbackward = tmpseries[~isnan(indpsdforward)]
rawpsd['indpsdbackward'] = indpsdbackward


#%%

### Save the created data to hdf5 file
savekw = {'complib':None,'complevel':0,'format':'fixed'}
## Surpressed for data conservation
# rawvar.to_hdf('hiwcdata.h5',key='rawvar',**savekw)
## save another rawvar as version hdf5 / netcdf4 using xarray.DataSet
# rawvar.index.name='time'
# xr.Dataset(rawvar).to_netcdf('pythondata/rawvards.h5',format='NETCDF4',mode='w')

# rawpsd.to_hdf('hiwcdata.h5',key='rawpsd',**savekw)
bin_div = np.arange(10.,12850.1,10.)
pd.DataFrame(bin_div).to_hdf('hiwcdata.h5',key='bin_div')
### rawpsd = pd.read_hdf('hiwcdata.h5',key='rawpsd') # For read in the file


#%%

### IGF routine
f_mygamma = lambda nml, x: 10**nml[0]*x**nml[1]*np.exp(-nml[2]*x)
def f_one_mode(psd, bin_div, moments):
    bin_diff = np.diff(bin_div)
    bin_mid = (bin_div[:-1]+bin_div[1:])/2.
    mobs = np.empty(moments.shape)
    for szmoment in range(len(mobs)):
        mobs[szmoment] = np.sum( psd*bin_diff*bin_mid**(moments[szmoment]) )

    x0 = np.array([log10(300), -1, 0.0014])
    tmpind = psd>0
    tmpfun = lambda x, nml1, nml2, nml3: f_mygamma(np.array([nml1,nml2,nml3]),x)
    nml0, pcov = scipy.optimize.curve_fit(tmpfun, bin_mid[tmpind], psd[tmpind], p0=x0)
    # nml0 = [1.8744, -0.6307, 0.0033]
    f_fit = lambda nml: f_sum_chisquare(nml, moments, bin_mid, bin_diff, mobs)
    optresult = scipy.optimize.minimize(f_fit, nml0, method='Nelder-Mead')
    return optresult
    
def f_sum_chisquare( nml, moments, bin_mid, bin_diff, mobs ):
    psd = f_mygamma(nml, bin_mid)
    ### We may drop this condition if the code works fine.
#     if any(np.isnan(psd)) or any(np.isinf(psd)):
#         return np.inf
    mfit = np.empty(mobs.shape)
    for szmoment in range(len(mfit)):
        mfit[szmoment] = np.sum( psd*bin_diff*bin_mid**(moments[szmoment]) )
    return np.sum( (mfit-mobs)**2/mobs/mfit )

## For median mass diameter calculation
def f_mmd(bin_div, msd):
    bin_diff = np.diff(bin_div)
    cmsd = np.concatenate( (np.array([0]),np.cumsum(msd*bin_diff)) )
    if cmsd[-1]<=0:
        return np.NaN
    cmsd /= cmsd[-1]
    indtmp = np.where(np.diff(cmsd>0.5)==1)[0]
    x1,x2,y1,y2 = bin_div[indtmp], bin_div[indtmp+1], cmsd[indtmp], cmsd[indtmp+1]
    mmd = (x2-x1)/(y2-y1)*(0.5-y1)+x1
    return mmd


#%%

## IGF fit and save the data
moments = np.array([0,2,3])
PSDs = rawpsd['psddmax'].as_matrix()
shp = PSDs.shape

output = [None]*shp[0]
for szi in range(shp[0]):
# for szi in range(6):
    psd = PSDs[szi,:]
    print('\r'+str(szi)+' ... ',end="")
    if np.sum(psd>0)>10:
        try:
            output[szi]=f_one_mode(psd,bin_div,moments)
            print('done',end="")
        except:
            output[szi]='Error'
            print('Error encountered')
    else:
        print('N/A',end="")

file = open('output.p', 'wb')
pickle.dump(output,file)
file.close()
print('Output file saved')


#%%

## Compute MMD
MSDs = rawpsd['msddmax'].as_matrix()
shp = MSDs.shape
MMD = np.empty(shp[0])
for szi in range(shp[0]):
# for szi in range(1):
#     msd = MSDs[20066,:]
    msd = MSDs[szi,:]
    MMD[szi] = f_mmd(bin_div, msd)

## Surpressed for data conservation
# file = open('mmd.p', 'wb')
# pickle.dump(MMD,file)
# file.close()


#%%

### Load the created data
rawvar = pd.read_hdf('pythondata/hiwcdata.h5',key='rawvar')
rawvards = xr.open_dataset('pythondata/rawvards.h5')
rawpsd = pd.read_hdf('pythondata/hiwcdata.h5',key='rawpsd')
bin_div = pd.read_hdf('pythondata/hiwcdata.h5',key='bin_div').as_matrix().ravel()

with open('pythondata/mmd.p', 'rb') as file:
    MMD = pickle.load(file)


#%%

## Transform psd data from pandas to xarray
bin_mid = (bin_div[1:]+bin_div[:-1])/2.
rawpsd.index.rename('time',inplace=True)
tmpda = []
for psdstr in rawpsd.keys().levels[0][:-1]:
    a=xr.DataArray(rawpsd[psdstr],coords=[('time',rawpsd.index),('bin',bin_mid)],name=psdstr)
    tmpda.append(a)
lastkey = rawpsd.keys().levels[0][-1]
a=xr.DataArray(rawpsd[lastkey],coords=[('time',rawpsd.index)],name=lastkey)
tmpda.append(a)
psdds = xr.merge(tmpda)

psdds.to_netcdf('pythondata/psdds.h5',format='NETCDF4',mode='w')


#%%

## combining rawvar with rasta and create large rastacombine.h5 file
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

print(datetime.datetime.now())
print('Reading rasta raw files ...')
tmplist = [loadRastaflt(szi) for szi in range(1,24)]
print(datetime.datetime.now())
print('Finished reading, combining ...')
# list(filter((None).__ne__, tmplist)) is a method to remove all the None cases in the list
rastadataset = xr.concat(list(filter((None).__ne__, tmplist)),dim='time')
print(datetime.datetime.now())
print('Writing to file ...')
rastadataset.to_netcdf('pythondata/rastacombine.h5',format='NETCDF4',mode='w')
print('Done.')
print(datetime.datetime.now())

