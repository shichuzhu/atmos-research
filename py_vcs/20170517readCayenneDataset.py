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
import netCDF4

import warnings
# warnings.filterwarnings('ignore')
# warnings.filterwarnings('default') # restore default settings

## For debug mode
from IPython.core.debugger import Tracer
#Tracer()() #this one triggers the debugger


#%%

# meta data for defining time range of flights
metatmp = [
    'F20_1Hz-HAIC-2015_core_v4_20150509_fs150009_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150510_fs150010_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150512_fs150011_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150514_fs150012_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150515_fs150013_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150516_fs150014_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150516_fs150015_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150518_fs150016_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150519_fs150017_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150523_fs150018_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150523_fs150019_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150524_fs150020_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150525_fs150021_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150526_fs150022_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150526_fs150023_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150527_fs150024_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150528_fs150025_20151106_V1.txt',
    'F20_1Hz-HAIC-2015_core_v4_20150529_fs150026_20151106_V1.txt']
    #                          26,    34,     41,43
# meta = [(tmp[26:34], int(tmp[41:43])) for tmp in metatmp]
meta = {int(tmp[41:43]):tmp[26:34] for tmp in metatmp}

outputfile = ''
# cayenne_bulk.nc  required for generating time-sync files
# cayenne_sync.nc  huge time sync file containing all data
# cayenne_sync_bulk.nc  all data except lamp psd
# cayenne_sync_gz.nc  zlib compressed lamp psd only


#%%

## For SAFFIRE V4 data
def read_saffire(filepath):
    tmp = os.path.basename(filepath)
    tmpflight = int(tmp[41:43])
    tmpdate = meta[tmpflight]
    '''
    Things to note
    1. The time is in seconds and can be non-sharp seconds.
    '''
    varnames = [
        'Timeinsecond',
        'event_marker',
        'latitude',
        'longitude',
        'altitude_gps',
        'altitude_airins',
        'platform_roll_angle',
        'platform_pitch_angle',
        'platform_orientation',
        'air_pressure',
        'air_temperature_rm',
        'air_temperature_impact',
        'air_temperature_adc',
        'dew_point_temperature',
        'relative_humidity',
        'humidity_mixing_ratio_hygrometer',
        'humidity_mixing_ratio_aerodata',
        'humidity_mixing_ratio_wvss2',
        'platform_speed_wrt_air',
        'platform_acceleration',
        'platform_course',
        'platform_speed_wrt_ground_aipov',
        'platform_course',
        'platform_speed_wrt_ground_gps',
        'upward_platform_speed_wrt_ground',
        'angle_of_attack',
        'angle_of_sideslip',
        'eastward_wind',
        'northward_wind',
        'upward_air_velocity',
        'wind_from_direction',
        'wind_speed',
        'mic_msofreqice_rs_sync_1'
    ]

    lookup = 'Warning : most measurements are not valid before take-off and after landing'
    comments = []
    with open(filepath) as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup in line:
                skipline = num
                break
    ### r"\s+" refers to one or more occurences of whitespace, while r"\s*" will match zero and would raise a warning
    tmp = pd.read_csv(filepath,skiprows=skipline,names=varnames,sep=r"\s+",na_values=3.40282347e+38)

    # pd.TimedeltaIndex(round(tmp.Timeinsecond).astype(int),units='s')
    tmp.loc[:,'time'] = (pd.TimedeltaIndex(tmp.Timeinsecond,unit='s') + pd.Timestamp(tmpdate)).round('s')
    tmp.set_index(keys='time',inplace=True,verify_integrity=False)
    tmp = tmp[pd.notnull(tmp.index)]
    tmp.index.rename('timeutc',inplace=True)

    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in saffire is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]
    tmp['flightnum'] = tmpflight
    return tmp

saffireset = [ read_saffire(filepath) for filepath 
          in glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataCayenne/20170513lamp/saffire/*') ]
# astype does not provide a 'coerce' option
datasetsaffire = pd.DataFrame().append(saffireset)
datasetsaffire.to_xarray().to_netcdf(outputfile,mode='a',
    format='NETCDF4',engine='netcdf4',group='/saffire',unlimited_dims=['timeutc'])
# datasetsaffire.apply(pd.to_numeric,errors='raise').dropna(how='all').equals(datasetsaffire) -> gives a True


#%%

# Read IKP
# with open('tmp/ikpdata.p','rb') as fin:
#     datasetipk = pickle.load(fin)

def read_ikp2(filepath):
    tmp = os.path.basename(filepath)
    tmpdate, tmpflight = tmp[4:8]+tmp[9:11]+tmp[12:14],int(tmp[23:25])
    # Skip the first row of each csv file. Warning this would lose data! Otherwise data are read as string, not float
    # tmp = pd.read_csv(filepath, header=4, skiprows=1, na_values=-999.0)
    tmp = pd.read_csv(filepath, header=4, na_values=-999.0)
    tmp.loc[:,'time'] = pd.Timestamp(tmpdate) + pd.to_timedelta(tmp['time'],errors='coerce') # Get rid of null lines
    tmp.set_index(keys='time',inplace=True,verify_integrity=False)
    tmp = tmp[pd.notnull(tmp.index)]
    tmp.index.rename('timeutc',inplace=True)

    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in IKP is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]
    return tmp

ikpset = [ read_ikp2(filepath) for filepath 
          in glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataCayenne/20170513lamp/ikp/*') ]

datasetipk = pd.DataFrame().append(ikpset).astype(float).replace(-999., nan)
datasetipk.to_xarray().to_netcdf(outputfile,mode='a', ## 'w' for first time, later change to 'a'
    format='NETCDF4',engine='netcdf4',group='/ikp',unlimited_dims=['timeutc'])


#%%

# Read robust
def read_robust(filepath):
    global meta
    tmp = os.path.basename(filepath)

    tmpflight = int(tmp[33:35])
    tmpdate = meta[tmpflight]
    # Skip the first row of each csv file. Warning this would lose data! Otherwise data are read as string, not float
    # tmp = pd.read_csv(filepath, header=4, skiprows=1, na_values=-999.0)
    tmp = pd.read_csv(filepath, header=0)
    tmp.columns = ['time','TWC_robust']
    tmp.loc[:,'time'] = pd.Timestamp(tmpdate) + pd.to_timedelta(tmp['time'],errors='coerce') # Get rid of null lines
    tmp.set_index(keys='time',inplace=True,verify_integrity=False)
    tmp = tmp[pd.notnull(tmp.index)]
    tmp.index.rename('timeutc',inplace=True)

    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in Robust is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]
    return tmp

robustset = [ read_robust(filepath) for filepath 
          in glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataCayenne/20170513lamp/robust/*V2*') ]
# astype does not provide a 'coerce' option
datasetrobust = pd.DataFrame().append(robustset).apply(pd.to_numeric, axis=0, errors='coerce').dropna()
datasetrobust.to_xarray().to_netcdf(outputfile,mode='a',
    format='NETCDF4',engine='netcdf4',group='/robust',unlimited_dims=['timeutc'])


#%%

## CDP Conc & PSD
def read_cdpconc(filepath):
    tmp = os.path.basename(filepath)
    tmpflight = int(tmp[24:26])
    tmpdate = meta[tmpflight]
    tmp = pd.read_csv(filepath)

    tmp.loc[:,'time'] = (pd.TimedeltaIndex(tmp.time,unit='s') + pd.Timestamp(tmpdate)).round('s')
    tmp.set_index(keys='time',inplace=True,verify_integrity=False)
    tmp = tmp[pd.notnull(tmp.index)]
    tmp.index.rename('timeutc',inplace=True)

    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in saffire is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]
    return tmp

tmpset = [ read_cdpconc(filepath) for filepath 
          in glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataCayenne/20170513lamp/cdp/Conc*') ]
# dropna applies changes
datasetcdpconc = pd.DataFrame().append(tmpset).apply(pd.to_numeric,errors='raise').dropna(how='all')
# datasetcdpconc.to_xarray().to_netcdf('tmp/cayenne_bulk.nc',mode='a',
#     format='NETCDF4',engine='netcdf4',group='/cdp',unlimited_dims=['timeutc'])

## CDP PSD
def read_cdppsd(filepath):
    tmp = os.path.basename(filepath)
    tmpflight = int(tmp[23:25])
    tmpdate = meta[tmpflight]
    tmp = pd.read_csv(filepath)

    tmp.loc[:,'time'] = (pd.TimedeltaIndex(tmp.time,unit='s') + pd.Timestamp(tmpdate)).round('s')
    tmp.set_index(keys='time',inplace=True,verify_integrity=False)
    tmp = tmp[pd.notnull(tmp.index)]
    tmp.index.rename('timeutc',inplace=True)

    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in saffire is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]
    return tmp

tmpset = [ read_cdppsd(filepath) for filepath 
          in glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataCayenne/20170513lamp/cdp/PSD*') ]
# dropna applies changes
datasetcdppsd = pd.DataFrame().append(tmpset).apply(pd.to_numeric,errors='raise').dropna(how='all')

# Combine and create xr.Dataset
b = xr.DataArray(datasetcdppsd)
b = b.rename({'dim_1':'bin_mid_cdp'})
b.bin_mid_cdp.values = b.bin_mid_cdp.astype(float)
b.bin_mid_cdp.attrs['unit']='um'
c= datasetcdpconc.to_xarray()
c['psd']=b
c.to_netcdf(outputfile,mode='a',
    format='NETCDF4',engine='netcdf4',group='/cdp',unlimited_dims=['timeutc'])


#%%

## LAMP PSD
def read_lamppsd(filepath,tmpflight,key):
    tmpdate = meta[tmpflight]
    tmp = pd.read_csv(filepath,sep=r"\s+",header=None,skiprows=1)

    tmp.iloc[:,0] = pd.TimedeltaIndex(tmp.iloc[:,0],unit='s') + pd.Timestamp(tmpdate)
    tmp.set_index(keys=0,inplace=True,verify_integrity=False)
    tmp = tmp[pd.notnull(tmp.index)]
    tmp.index.rename('timeutc',inplace=True)

    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in lamppsd is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]
    
    tmp = pd.concat([tmp], axis=1, keys=[key]) # Add the outmost row index flight number

    return tmp

# Generate fileinfo for lamppsd
lamproot = '/data/mcfarq/a/szhu28/research/HIWC/data/fulldataCayenne/20170513lamp/dataHAIC/'
paras = np.array([
    ['msddeq','msddmax','psddeq','psddmax','psdly'],
    ['*MassSize*Deq*','*MassSize*Dmax*','*Composite*Deq*','*Composite*Dmax*','*Composite*Ly*']])

tmpinfo = []
for i in range(paras.shape[1]):
    wildcard = paras[1][i]
    filepaths = glob.glob(lamproot+wildcard)
    filepaths.sort()
    tmpset = []
    for filepath in filepaths:
        tmp = os.path.basename(filepath)
        m = re.search('(?<=_VOL)\d\d(?=-)', tmp)
        tmpflight = int(m.group(0))
        tmpinfo.append({'filepath':filepath,'tmpflight':tmpflight,'key':paras[0][i]})
fileinfo = np.array(tmpinfo).reshape(paras.shape[1],-1)

# In the sense of psd binned
tmpdataset = []
for i in range(fileinfo.shape[0]):
    tmpset = []
    for j in range(fileinfo.shape[1]):
        tmpset.append(read_lamppsd(**fileinfo[i][j]))
    tmp = pd.DataFrame().append(tmpset)
    tmpdataset.append( tmp )

b = pd.concat(tmpdataset,axis=1)
b = pd.Panel({x:b[x] for x in b.columns.levels[0]}).to_xarray()

b = b.rename({'major_axis':'timeutc','minor_axis':'bin_mid_composite'})
b = b.to_dataset(dim='items')
b.to_netcdf(outputfile,mode='a',
    format='NETCDF4',engine='netcdf4',group='/lamp',unlimited_dims=['timeutc'])


#%%

# Merge all data to one time dimension
setnamelist = ['saffire','ikp','robust','cdp','lamp']
# setnamelist = ['saffire','ikp','robust','cdp']
alldata = []
for name in setnamelist:
    alldata.append( xr.open_dataset('tmp/cayenne_bulk.nc',group=name) )

alltimelist = [x.timeutc for x in alldata]
alltime = unique(xr.concat(alltimelist,dim='timeutc'))
alltime = xr.DataArray(alltime,dims={'timeutc':len(alltime)},coords={'timeutc':alltime})
alltime = alltime.to_dataset(name='time_to_drop')

outputfile = 'tmp/cayenne_sync_bulk.nc'
for i in range(len(setnamelist)-1):
    name = setnamelist[i]
#     mode = 'a'
    if i == 0:
        mode = 'w'
    else:
        mode = 'a'
    xr.merge([alltime,alldata[i]]).drop('time_to_drop').to_netcdf(outputfile,mode=mode,
    format='NETCDF4',engine='netcdf4',group='/'+name,unlimited_dims=['timeutc'])

i = 4
name = setnamelist[i]
outputfile = 'tmp/cayenne_sync_gz.nc'
xr.merge([alltime,alldata[i]]).drop('time_to_drop').to_netcdf(outputfile,mode='w',
    format='NETCDF4',engine='netcdf4',group='/'+name,unlimited_dims=['timeutc'],
    encoding = {x:{'zlib':True} for x in ['msddeq','msddmax','psddeq','psddmax','psdly']} )


#%%

# Elaborate coordinates and attributes
totalbins = 1284
bin_div = 10*np.arange(totalbins+1)+10
bin_mid = (bin_div[:-1]+bin_div[1:])/2
with netCDF4.Dataset('tmp/cayenne_sync_gz.nc','a') as dset:
    dset['/lamp/bin_mid_composite'][:] = bin_mid
    dset['/lamp'].bin_div = bin_div.astype(float)
    dset['/lamp'].bin_div_units = 'um'


#%%

### IGF routine
f_mygamma = lambda nml, x: 10**nml[0]*x**nml[1]*np.exp(-nml[2]*x)
def f_one_mode(psd, bin_div, moments):
    # bin_diff = np.diff(bin_div)
    # bin_mid = (bin_div[:-1]+bin_div[1:])/2.
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
    # bin_diff = np.diff(bin_div)
    cmsd = np.concatenate( (np.array([0]),np.cumsum(msd*bin_diff)) )
    if cmsd[-1]<=0:
        return np.NaN
    cmsd /= cmsd[-1]
    indtmp = np.where(np.diff(cmsd>0.5)==1)[0]
    x1,x2,y1,y2 = bin_div[indtmp], bin_div[indtmp+1], cmsd[indtmp], cmsd[indtmp+1]
    mmd = (x2-x1)/(y2-y1)*(0.5-y1)+x1
    return mmd


#%%

# IGF fit
a = xr.open_dataset('tmp/cayenne_sync_gz.nc',group='/lamp')
bin_div = a.bin_div
bin_diff = np.diff(bin_div)
bin_mid = (bin_div[:-1]+bin_div[1:])/2.

## IGF fit and save the data
moments = np.array([0,2,3])
PSDs = a['psddmax'].values
MSDs = a['msddmax'].values
shp = PSDs.shape
# At least 11 non-zero bins are required for a valid PSD
validpsdbool = (~any(isnan(PSDs),axis=1)) & (sum(PSDs>0,axis=1)>10)
validpsdind = np.where(validpsdbool)[0]

'''
import signal
class TimeoutException(Exception):   # Custom exception class
    pass
def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException
# Change the behavior of SIGALRM
signal.signal(signal.SIGALRM, timeout_handler)

for i in range(3):
    # Start the timer. Once 5 seconds are over, a SIGALRM signal is sent.
    signal.alarm(5)    
    # This try/except loop ensures that 
    #   you'll catch TimeoutException when it's sent.
    try:
        A(i) # Whatever your function that might hang
    except TimeoutException:
        continue # continue the for loop if function A takes more than 5 second
    else:
        # Reset the alarm
        signal.alarm(0)
'''
# IGF fit and save data

output = {}
i = 0
totali = len(validpsdind)
for szi in validpsdind:
    i+=1
    psd = PSDs[szi,:]
    print(str(szi)+' ... '+str(i/totali),end="")
    try:
        output[szi]=f_one_mode(psd,bin_div,moments)
        print('\r'+str(szi)+' Done. Now',end="")
    except:
        output[szi]=None
        print('\r'+str(szi)+' Error. Now',end="")

with open('tmp/output.p', 'wb') as file:
    pickle.dump(output,file)
    print('Output file saved')


#%%

# Read the output fit file and save to cayenne_sync_bulk.nc as /lampproc group
with open('tmp/output.p', 'rb') as file:
    output = pickle.load(file)
    print('Output file read')

with netCDF4.Dataset("tmp/cayenne_sync_gz.nc", "r", format="NETCDF4") as file:
    varname = list(file.groups['lamp'].variables.keys())
    varname.remove('timeutc')
frame = xr.open_dataset('tmp/cayenne_sync_gz.nc',group='/lamp',drop_variables=varname)

tmpnml = np.empty((frame.dims['timeutc'],3), dtype=float)
tmpnml[:] = nan
output = {k:v.x for k,v in output.items() if v is not None}
tmpnml[ list(output.keys()),:] = np.array([x for x in output.values()])
tmpvalidbinnum = np.empty(PSDs.shape[0])
tmpvalidbinnum[:] = nan
mask = ~any(isnan(PSDs),axis=1)
tmpvalidbinnum[mask] = np.sum(PSDs[mask,:]>0,axis=1)
frame['validbinnum'] = xr.DataArray(tmpvalidbinnum,dims=['timeutc'])
frame['nml'] = xr.DataArray(tmpnml,dims=['timeutc','dimnml'])

frame.to_netcdf('tmp/cayenne_sync_bulk.nc',mode='a',
    format='NETCDF4',engine='netcdf4',group='/lampproc',unlimited_dims=['timeutc'])


#%%

# Add mmd to cayenne_sync_bulk.nc, also as an template for adding variable
tmpmmd = np.empty(PSDs.shape[0])
tmpmmd[:] = nan
validMSDs = MSDs[validpsdind,:]
validmmd = np.array([f_mmd(bin_div,x) for x in validMSDs])
tmpmmd[validpsdind] = validmmd

# As a template to add new variable to existing netcdf file instead of using xarray
file = netCDF4.Dataset('tmp/cayenne_sync_bulk.nc',mode='a')
grp = file['/lampproc']
# The extra comma makes sure the passed constant is a tuple as required
# Once variable created, the file will reject repeated creation
# varmmd = grp.createVariable("mmd","f8",('timeutc',))
varmmd = grp['mmd']
varmmd[:] = tmpmmd
file.close()


#%%

# Destructively concatenate RASTA data to the merged dataset
outputfile = 'tmp/cayenne_rasta.nc'
# Read RASTA files
def read_rasta(filepath):
    filename = os.path.basename(filepath)
    tmpflight = int(re.search('(?<=_F)\d*(?=_radonvar)',filename).group())
    rasta = xr.open_dataset(filepath)

    sec = 3600*rasta.time.values
    time = (pd.TimedeltaIndex(sec,unit='s') + pd.Timestamp(meta[tmpflight])).round('s')
    rasta['time'] = time

    if time.is_unique is False:
        print('Warning!! Duplicate data in saffire is found, dropping ->')
        print(sum(time.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        uniindex = np.where(~time.duplicated(keep='last'))[0]
        rasta = rasta.isel(time=uniindex)
    
    rasta.rename({'time':'timeutc'},inplace=True)
    return rasta

rastafilenames = glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataCayenne/20170513lamp/rasta/*.nc')

rastaset = [ read_rasta(filepath) for filepath in rastafilenames ]
timeutc = xr.open_dataset('tmp/cayenne/cayenne_sync_bulk.nc',group='/saffire').timeutc

key = 'timeutc'
xr.Dataset(data_vars={'timeutc':timeutc}).to_netcdf(outputfile,mode='w', ## 'w' for first time, later change to 'a'
    format='NETCDF4',engine='netcdf4',group='/',unlimited_dims=['timeutc'])

keys = list(rastaset[0].data_vars.keys())
tot = len(keys)
i = 0
for key in keys:
    i+=1
    print('Processing key '+key+' ... '+str(i)+'/'+str(tot))
    tmpset = [ x[key] for x in rastaset ]

    datasettmp = xr.concat(tmpset,dim='timeutc')
    datasettmp = xr.align(datasettmp,indexes={'timeutc':timeutc})[0]
    datasettmp.to_netcdf(outputfile,mode='a', ## 'w' for first time, later change to 'a'
        format='NETCDF4',engine='netcdf4',group='/'+key,unlimited_dims=['timeutc'],
        encoding={ key:{'zlib':True, 'complevel':1} } )


#%%

# Construct link table to MTSAT file
dirpaths = glob.glob('/data/gpm/a/shared/szhu28/hiwcproc/cayenne/mtsat/*')
dirpaths.sort()

allfiles = []
for dirpath in dirpaths:
    datestr = os.path.basename(dirpath)
    filepaths = glob.glob(dirpath+'/*')
    filepaths.sort()
    for filepath in filepaths:
        filename = os.path.basename(filepath)
        timestr = os.path.splitext(filename)[0]
        allfiles.append((datestr+timestr,filepath))

timesat = pd.DatetimeIndex( [ x[0] for x in allfiles ] )
satfiles = np.array( [ x[1] for x in allfiles ] )
satfiles = xr.DataArray(satfiles, dims=['timesat'], coords={'timesat':timesat}, name='satfiles')

sf = xr.open_dataset('tmp/cayenne/cayenne_sync_bulk.nc',group='/saffire')[ ['latitude','longitude'] ]
timeutc = sf.timeutc.values
orig = satfiles.timesat.values
insind = np.searchsorted(orig,timeutc)
torf = abs(timeutc-orig[insind-1])>abs(timeutc-orig[insind])
inserted = np.where(torf, orig[insind], orig[insind-1])

indarr = np.where(np.insert(diff(inserted),0,0))[0]
indarr = np.concatenate([[0],indarr,[None]])

varkeys = [
 'latitude',
 'longitude',
 'reflectance_vis',
 'visible_count',
 'reflectance_nir',
 'temperature_sir',
 'temperature_67',
 'temperature_ir',
 'temperature_sw',
 'broadband_shortwave_albedo',
 'broadband_longwave_flux',
 'cloud_ir_emittance',
 'cloud_phase',
 'cloud_visible_optical_depth',
 'cloud_particle_size',
 'cloud_lwp_iwp',
 'cloud_effective_temperature',
 'cloud_top_pressure',
 'cloud_effective_pressure',
 'cloud_bottom_pressure',
 'cloud_top_height',
 'cloud_effective_height',
 'cloud_bottom_height',
 'cloud_top_temperature',
 'cloud_bottom_temperature',
 'pixel_skin_temperature',
]

satds = xr.Dataset(data_vars={'timeutc':timeutc})

def f_dist(x1,y1,x2,y2):
    # remember to switch to radial before using
    x1,y1,x2,y2 = [ f_rad(x) for x in [x1,y1,x2,y2] ]
    return np.sqrt(((x1-x2)*np.cos((y1+y2)/2))**2 + (y1-y2)**2)
def f_rad(x):
    return x/180*np.pi

def f_4ptinterp(xi,yi,xs,ys,zs):
    dists = f_dist(xi,yi,xs,ys)
    return np.average(zs,weights=np.minimum(1/dists,1e10))


#%%

# Process each MTSAT files and generate along-flight dataset
toconcat = []
outputfile = 'tmp/cayenne_mtsat.nc'
warnings.filterwarnings('ignore')
for i in range(len(indarr)-1):
    print(i)
# for i in range(2):
    indi = slice(indarr[i],indarr[i+1])
    filepath = str(satfiles.sel(timesat=inserted[indarr[i]]).values)

    # Gladly in Cayenne dataset missing_value attribute works fine compared to Darwin
    tmpds = xr.open_dataset(filepath, mask_and_scale=True)
    tmpds = tmpds[varkeys]

    frameds = sf.isel(timeutc=indi)

    lons = tmpds['longitude'].values
    lats = tmpds['latitude'].values
    lon = frameds.longitude.values
    lat = frameds.latitude.values

    # Generate two arrays for lon, lat
    # 1. smallest value greater than all the elements before (inclusive) the current line (sgtb)
    # 2. greatest value smaller than all the elements after (inclusive) the current line (gsta)

    # For latitude, the generate trend is DECREASING wrt. index increasing.
    tmp = np.nanmax(lons,axis=0)
    tmp[isnan(tmp)] = -inf
    tmp = np.maximum.accumulate(tmp)
    tmpind = np.nonzero(tmp==-inf)[0][-1]
    tmp[:tmpind+1] = tmp[tmpind+1]
    sgtblon = tmp

    tmp = np.nanmin(lons,axis=0)
    tmp[isnan(tmp)] = inf
    tmp = np.minimum.accumulate(tmp[::-1])[::-1]
    tmpind = np.nonzero(tmp==inf)[0][0]
    tmp[tmpind:] = tmp[tmpind-1]
    gstalon = tmp

    # any(gstalon>sgtblon) This should be false if the above codes work fine.

    # For latitude, the generate trend is DECREASING wrt. index increasing.
    tmp = np.nanmax(lats,axis=1)
    tmp[isnan(tmp)] = -inf
    tmp = np.maximum.accumulate(tmp[::-1])[::-1]
    tmpind = np.nonzero(tmp==-inf)[0][0]
    tmp[tmpind:] = tmp[tmpind-1]
    sgtblat = tmp

    tmp = np.nanmin(lats,axis=1)
    tmp[isnan(tmp)] = inf
    tmp = np.minimum.accumulate(tmp)
    tmpind = np.nonzero(tmp==inf)[0][-1]
    tmp[:tmpind+1] = tmp[tmpind+1]
    gstalat = tmp
    # any(gstalat>sgtblat) This should be false if the above codes work fine.

    # gsta -> upper bound, sgtb -> lower bound
    lonr = np.searchsorted(gstalon, lon)+1 # +1 is for the upper bound exclusive in python
    lonl = np.searchsorted(sgtblon, lon)-1 # -1 is for considering all possibility

    # gsta -> upper bound, sgtb -> lower bound
    latl = len(gstalat)-np.searchsorted(gstalat[::-1], lat)-1 # -1 is for considering all possibility
    latr = len(sgtblat)-np.searchsorted(sgtblat[::-1], lat)+1 # +1 is for the upper bound exclusive in python

    # Add another wrapper for key values loop
    toaddds = {}
    for key in varkeys:
        tmpzs = tmpds[key].values
        zint = []
        for j in range(len(lon)):
            indsubgrid = slice(latl[j],latr[j]),slice(lonl[j],lonr[j])
            x = lons[indsubgrid]
            if len(x) == 0:
                zint.append(nan)
                continue
            y = lats[indsubgrid]
            z = tmpzs[indsubgrid]
            indmin = np.unravel_index( ((x-lon[j])**2+(y-lat[j])**2).argmin(), x.shape)

            indnine = slice(indmin[0]-1,indmin[0]+2),slice(indmin[1]-1,indmin[1]+2)
            xs = x[indnine].ravel()
            if len(xs) == 0:
                zint.append(nan)
                continue
            ys = y[indnine].ravel()
            zs = z[indnine].ravel()
            zint.append(f_4ptinterp(lon[j],lat[j],xs,ys,zs))
        zint = np.array(zint)
        toaddds[key] = xr.DataArray(zint, dims=['timeutc'], coords={'timeutc':frameds.timeutc},
                                    attrs=tmpds[key].attrs, name=key)
    toconcat.append(xr.Dataset(toaddds))

warnings.filterwarnings('default') # restore default settings

mtsatproc = xr.concat(toconcat,dim='timeutc')

# Remember to add the variable for time difference and a mask
tmpmtsatproc = xr.Dataset({},coords={'timeutc':timeutc})
tmp = timeutc - inserted
tmpmtsatproc['timelag'] = xr.DataArray(tmp, dims=['timeutc'], coords={'timeutc':timeutc})
# Note the threshold is 30 min for Darwin and 15 min for Cayenne
tmpmtsatproc['validlagmask'] = xr.DataArray( tmp <= np.timedelta64(15,'m') , dims=['timeutc'], coords={'timeutc':timeutc})
mtsatproc = xr.merge([tmpmtsatproc,mtsatproc],join='left')
mtsatproc.to_netcdf(outputfile,mode='w',
    format='NETCDF4',engine='netcdf4',group='/mtsatproc',unlimited_dims=['timeutc'])


#%%

# Convert and create NCAR ECMWF archive
import pygrib
import resource

# Parameters for Cayenne
filepaths = glob.glob('/data/gpm/a/shared/szhu28/hiwcproc/era_interim/2015**/*regn128sc*',recursive=True)
kw = {'lat1':0,'lat2':10,'lon1':-60+360,'lon2':-45+360}
outputnoext = 'tmp/cayenne/ecmwf'

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

