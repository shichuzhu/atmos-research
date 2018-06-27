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

# Check memory usage
import resource
print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)


#%%

# meta data for defining time range of flights
metatmp = [
    #                           27,    35,     42,44
    'F20_1Hz-HAIC_base_aipov_v4_20140116_fs140001.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140116_fs140002.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140117_fs140003.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140118_fs140004.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140121_fs140005.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140123_fs140006.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140124_fs140007.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140127_fs140008.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140128_fs140009.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140129_fs140010.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140130_fs140011.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140202_fs140012.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140203_fs140013.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140204_fs140014.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140205_fs140015.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140207_fs140016.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140208_fs140017.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140208_fs140018.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140209_fs140019.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140210_fs140020.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140217_fs140021.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140217_fs140022.txt',
    'F20_1Hz-HAIC_base_aipov_v4_20140218_fs140023.txt']

# meta = [(tmp[26:34], int(tmp[41:43])) for tmp in metatmp]
meta = {int(tmp[42:44]):tmp[27:35] for tmp in metatmp}

outputfile = 'tmp/darwin_bulk.nc'
# cayenne_bulk.nc  required for generating time-sync files
# cayenne_sync.nc  huge time sync file containing all data
# cayenne_sync_bulk.nc  all data except lamp psd
# cayenne_sync_gz.nc  zlib compressed lamp psd only


#%%

## For SAFFIRE V4 data
def read_saffire(filepath):
    tmp = os.path.basename(filepath)
    tmpflight = int(tmp[42:44])
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
filelist = glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/data_src/data20161012/saffire/F20_1Hz-*')
filelist.sort()
saffireset = [ read_saffire(filepath) for filepath in filelist ]
# astype does not provide a 'coerce' option
datasetsaffire = pd.DataFrame().append(saffireset)
datasetsaffire.to_xarray().to_netcdf(outputfile,mode='w',
    format='NETCDF4',engine='netcdf4',group='/saffire',unlimited_dims=['timeutc'])
# datasetsaffire.apply(pd.to_numeric,errors='raise').dropna(how='all').equals(datasetsaffire) -> gives a True

# Create metaspan useful for files not named by flight number
test = datasetsaffire.flightnum
def tmpfun(flt):
    tmp = test[test == flt].index
    return [flt,min(tmp),max(tmp)]
metaspan = list(map(tmpfun,range(1,24)))
takeofftime = pd.DatetimeIndex([x[1] for x in metaspan])


#%%

# Read IKP
# with open('tmp/ikpdata.p','rb') as fin:
#     datasetipk = pickle.load(fin)

def read_ikp2(filepath):
    tmp = os.path.basename(filepath)
    ts = tmp
    tsflt = pd.Timestamp(ts[4:8]+ts[9:11]+ts[12:14]+ts[15:21])
    tmpflight = argmin(abs(tsflt - takeofftime))+1
        
    tmpdate = meta[tmpflight]
    # Skip the first row of each csv file. Warning this would lose data! Otherwise data are read as string, not float
    # tmp = pd.read_csv(filepath, header=4, skiprows=1, na_values=-999.0)
    
    # Weird line 7 ending comma bug leads to discarding first two lines of data.
    delrow = list(range(7))
    delrow.remove(4)
    tmp = pd.read_csv(filepath, header=0, skiprows=delrow, na_values=-999.0)
    
    tmp.loc[:,'Stimech'] = pd.Timestamp(tmpdate) + pd.to_timedelta(tmp['Stimech'],errors='coerce') # Get rid of null lines
    tmp.set_index(keys='Stimech',inplace=True,verify_integrity=False)
    tmp = tmp[pd.notnull(tmp.index)]
    tmp.index.rename('timeutc',inplace=True)

    if tmp.index.is_unique is False:
        print('Warning!! Duplicate data in IKP is found, dropping ->')
        print(sum(tmp.index.duplicated()))
        ### Dropping duplicate based on index, see
        ### http://pandas.pydata.org/pandas-docs/stable/indexing.html#duplicate-data
        tmp = tmp[~tmp.index.duplicated(keep='last')]
    return tmp
filelist = glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/data_src/data20161012/ikp/f20-2014*')
filelist.sort()
ikpset = [ read_ikp2(filepath) for filepath in filelist ]

datasetipk = pd.DataFrame().append(ikpset).astype(float).replace(-999., nan)
datasetipk.to_xarray().to_netcdf(outputfile,mode='a', ## 'w' for first time, later change to 'a'
    format='NETCDF4',engine='netcdf4',group='/ikp',unlimited_dims=['timeutc'])


#%%

# Read robust
def read_robust(filepath):
    global meta
    tmp = os.path.basename(filepath)

    tmpflight = int(tmp[15:17])
    tmpdate = meta[tmpflight]
    # Skip the first row of each csv file. Warning this would lose data! Otherwise data are read as string, not float
    # tmp = pd.read_csv(filepath, header=4, skiprows=1, na_values=-999.0)
    tmp = pd.read_excel(filepath)
    tmp.columns = ['time','TWC_robust']
    # excel files have date and time ready
    tmp.loc[:,'time'] = pd.DatetimeIndex(tmp['time'],errors='coerce').round('1S') # Get rid of null lines
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

filelist = glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/dataHAIC/Robust_data_flt*')
filelist.sort()
robustset = [ read_robust(filepath) for filepath in filelist ]
# astype does not provide a 'coerce' option
datasetrobust = pd.DataFrame().append(robustset).apply(pd.to_numeric, axis=0, errors='coerce').dropna()
datasetrobust.to_xarray().to_netcdf(outputfile,mode='a',
    format='NETCDF4',engine='netcdf4',group='/robust',unlimited_dims=['timeutc'])


#%%

## CDP Conc
def read_cdpconc(filepath,tmpflight):
    tmp = os.path.basename(filepath)
    tmpdate = meta[tmpflight]
    tmp = pd.read_excel(filepath)
    tmp.columns = ['time','cdpconcpercm3']

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

# Need regular expression for flight number as they are not fixed width
tmpfileset = [ ( int(re.search('(?<=CDP3V)\d*(?=\.xls)', os.path.basename(filepath)).group()), filepath ) for filepath 
          in glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/dataHAIC/HAIC2014-Data-CDP3*') ]
tmpfileset.sort()
tmpset = [ read_cdpconc(fileinfo[1],fileinfo[0]) for fileinfo in tmpfileset ]
# dropna applies changes
datasetcdpconc = pd.DataFrame().append(tmpset).apply(pd.to_numeric,errors='raise').dropna(how='all')
datasetcdpconc.to_xarray().to_netcdf(outputfile,mode='a',
    format='NETCDF4',engine='netcdf4',group='/cdp',unlimited_dims=['timeutc'])


#%%

# Generate sorted fileinfo for lamp PSD & MSD. keyword: re regular expression
def tmpfun(filepath):
    filename = os.path.basename(filepath)
    tmpflight = int(re.search('(?<=vol)\d*(?=\.txt)',filename).group())
    # The (?:.(?!-))* is for non-greedy matching, reference http://stackoverflow.com/a/2527791/5426033
    signit = re.search('(?<=-)(?:.(?!-))*(?=_vol)',filename).group()
    key = dict(zip(['Deq','Dmax','ly'],['psddeq','psddmax','psdly']))[signit]
    return (key,tmpflight,filepath)
finfo1 = list(map(tmpfun,glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/lamp/Composite_2DS-PIP*')))

def tmpfun(filepath):
    filename = os.path.basename(filepath)
    tmpflight = int(re.search('(?<=F#)\d*(?=-Mass)',filename).group())
    signit = re.search('(?<=SizeD_).*(?=_2DS)',filename).group()
    key = dict(zip(['Deq','Dmax'],['msddeq','msddmax']))[signit]
    return (key,tmpflight,filepath)
finfo2 = list(map(tmpfun,glob.glob(
    '/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/data_src/data20151203/dataHAIC/*MassSize*')))
finfo = finfo1 + finfo2
finfo.sort()


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

bunch = [ read_lamppsd(x[2],x[1],x[0]) for x in finfo ]

lamppd = pd.concat( [ pd.DataFrame().append([bunch[finfo.index(x)] for x in finfo if x[0] == y ]) 
       for y in set([x[0] for x in finfo]) ] , axis=1 )

# with open('tmp/lamptmp.p','wb') as fout:
#     pickle.dump(lamppd,fout)

# with open('tmp/lamptmp.p','rb') as fin:
#     lamppd = pickle.load(fin)

# Since ly has different bin_div against other PSD, it is taken out to deal with alone.
partlist = list(lamppd.columns.levels[0])
partlist.remove('psdly')

b = pd.Panel({x:lamppd[x] for x in partlist}).to_xarray()
b = b.rename({'major_axis':'timeutc','minor_axis':'bin_mid_composite'}).to_dataset(dim='items')

b['bin_mid_composite'] = (1+np.arange(len(b.bin_mid_composite)))*10+5.
b.attrs['bin_div_composite'] = (1+np.arange(1+len(b.bin_mid_composite)))*10.
bt = b

b = pd.Panel({'psdly':lamppd['psdly']}).to_xarray()
b = b.rename({'major_axis':'timeutc','minor_axis':'bin_mid_ly'}).to_dataset(dim='items')

fn = '/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/lamp/Composite_2DS-PIP-Intercomp_PSD-ly_vol01.txt'
titletmp = pd.read_table(fn)
b['bin_mid_ly'] = np.array(titletmp.columns.values)[1:].astype(float)
bin_mid = np.array(b['bin_mid_ly'])
bin_div = np.append(np.empty_like(bin_mid),nan)
bin_div[0] = 10
for i in range(len(bin_mid)):
    bin_div[i+1] = 2*bin_mid[i]-bin_div[i]
b.attrs['bin_div_ly'] = bin_div

datasetlamp = xr.merge([b,bt])
datasetlamp.to_netcdf(outputfile,mode='a',
    format='NETCDF4',engine='netcdf4',group='/lamp',unlimited_dims=['timeutc'])


#%%

# Merge all data to one time dimension
outputfile = 'tmp/darwin_sync_bulk.nc'
outputfile1 = 'tmp/darwin_sync_gz.nc'
setnamelist = ['saffire','ikp','robust','cdp','lamp']
# setnamelist = ['saffire','ikp','robust','cdp']
alldata = []
for name in setnamelist:
    alldata.append( xr.open_dataset('tmp/darwin_bulk.nc',group=name) )

alltimelist = [x.timeutc for x in alldata]
alltime = unique(xr.concat(alltimelist,dim='timeutc'))
alltime = xr.DataArray(alltime,dims={'timeutc':len(alltime)},coords={'timeutc':alltime})
alltime = alltime.to_dataset(name='time_to_drop')

for i in range(len(setnamelist)-1):
    name = setnamelist[i]
#     mode = 'a'
    if i == 0:
        mode = 'w'
    else:
        mode = 'a'
    mergebulk = xr.merge([alltime,alldata[i]]).drop('time_to_drop')
    mergebulk.to_netcdf(outputfile,mode=mode,
    format='NETCDF4',engine='netcdf4',group='/'+name,unlimited_dims=['timeutc'])

i = 4
name = setnamelist[i]
mergepsd = xr.merge([alltime,alldata[i]]).drop('time_to_drop')
mergepsd.to_netcdf(outputfile1,mode='w',
    format='NETCDF4',engine='netcdf4',group='/'+name,unlimited_dims=['timeutc'],
    encoding = {x:{'zlib':True} for x in ['msddeq','msddmax','psddeq','psddmax','psdly']} )


#%%

# Elaborate coordinates and attributes
totalbins = 1284
bin_div = 10*np.arange(totalbins+1)+10
bin_mid = (bin_div[:-1]+bin_div[1:])/2

# outputfile1 = 'tmp/darwin_sync_gz.nc' created in the previous cell
with netCDF4.Dataset(outputfile1,'a') as dset:
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
        mobs[szmoment] = np.nansum( psd*bin_diff*bin_mid**(moments[szmoment]) )

    x0 = np.array([log10(300), -1, 0.0014])
    # any zero bins as well as nan will be ignored
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
#     if len(indtmp)>1:
#         print(indtmp)
    x1,x2,y1,y2 = bin_div[indtmp], bin_div[indtmp+1], cmsd[indtmp], cmsd[indtmp+1]
    mmd = (x2-x1)/(y2-y1)*(0.5-y1)+x1
    return mmd


#%%

get_ipython().run_cell_magic('capture', '--no-stdout', '# IGF fit\n# outputfile = \'tmp/darwin_sync_bulk.nc\'\n# outputfile1 = \'tmp/darwin_sync_gz.nc\'\na = xr.open_dataset(outputfile1,group=\'/lamp\')\nbin_div = a.bin_div\nbin_diff = np.diff(bin_div)\nbin_mid = (bin_div[:-1]+bin_div[1:])/2.\n\n## IGF fit and save the data\nmoments = np.array([0,2,3])\nPSDs = a[\'psddmax\'].values\nMSDs = a[\'msddmax\'].values\nshp = PSDs.shape\n# At least 11 non-zero bins are required for a valid PSD\n# Since last few bins of PSD is alwasys nan, use \'not all are nan\' as nalid PSD index\nvalidpsdbool = (~all(isnan(PSDs),axis=1)) & (sum(PSDs>0,axis=1)>10)\nvalidpsdind = np.where(validpsdbool)[0]\n\n\'\'\'\nimport signal\nclass TimeoutException(Exception):   # Custom exception class\n    pass\ndef timeout_handler(signum, frame):   # Custom signal handler\n    raise TimeoutException\n# Change the behavior of SIGALRM\nsignal.signal(signal.SIGALRM, timeout_handler)\n\nfor i in range(3):\n    # Start the timer. Once 5 seconds are over, a SIGALRM signal is sent.\n    signal.alarm(5)    \n    # This try/except loop ensures that \n    #   you\'ll catch TimeoutException when it\'s sent.\n    try:\n        A(i) # Whatever your function that might hang\n    except TimeoutException:\n        continue # continue the for loop if function A takes more than 5 second\n    else:\n        # Reset the alarm\n        signal.alarm(0)\n\'\'\'\n# IGF fit and save data\n\noutput = {}\ni = 0\ntotali = len(validpsdind)\nfor szi in validpsdind:\n    i+=1\n    psd = PSDs[szi,:]\n    print(str(szi)+\' ... \', i, totali, \'{:.2f}\'.format(i/totali*100),end="")\n    try:\n        output[szi]=f_one_mode(psd,bin_div,moments)\n        print(\'\\r\'+str(szi)+\' Done. Now\',end="")\n    except:\n        output[szi]=None\n        print(\'\\r\'+str(szi)+\' Error. Now\',end="")\nprint(\'Finished!\')\nwith open(\'tmp/output.p\', \'wb\') as file:\n    pickle.dump(output,file)\n    print(\'Output file saved\')')


#%%

# Read the output fit file and save to cayenne_sync_bulk.nc as /lampproc group
outputfile = 'tmp/darwin_sync_bulk.nc'
outputfile1 = 'tmp/darwin_sync_gz.nc'
with open('tmp/output.p', 'rb') as file:
    output = pickle.load(file)
    print('Output file read')

with netCDF4.Dataset(outputfile1, "r", format="NETCDF4") as file:
    varname = list(file.groups['lamp'].variables.keys())
    varname.remove('timeutc')
frame = xr.open_dataset(outputfile1,group='/lamp',drop_variables=varname)

tmpnml = np.empty((frame.dims['timeutc'],3), dtype=float)
tmpnml[:] = nan
# Remove entries that are None
output = {k:v.x for k,v in output.items() if v is not None}
tmpnml[ list(output.keys()),:] = np.array([x for x in output.values()])
tmpvalidbinnum = np.empty(PSDs.shape[0])
tmpvalidbinnum[:] = nan
mask = ~all(isnan(PSDs),axis=1)
tmpvalidbinnum[mask] = np.sum(PSDs[mask,:]>0,axis=1)
frame['validbinnum'] = xr.DataArray(tmpvalidbinnum,dims=['timeutc'])
frame['nml'] = xr.DataArray(tmpnml,dims=['timeutc','dimnml'])

frame.to_netcdf(outputfile,mode='a',
    format='NETCDF4',engine='netcdf4',group='/lampproc',unlimited_dims=['timeutc'])


#%%

# Add mmd to cayenne_sync_bulk.nc, also as an template for adding variable

# find valid msd index here. MSD is not exactly same as PSD in Darwin dataset
validmsdbool = (~all(isnan(MSDs),axis=1)) & (sum(MSDs>0,axis=1)>10)
validmsdind = np.where(validmsdbool)[0]

tmpmmd = np.empty(PSDs.shape[0])
tmpmmd[:] = nan
validMSDs = MSDs[validmsdind,:]
validmmd = np.array([f_mmd(bin_div,x) for x in validMSDs])
tmpmmd[validmsdind] = validmmd

# As a template to add new variable to existing netcdf file instead of using xarray
file = netCDF4.Dataset(outputfile,mode='a')
grp = file['/lampproc']
# The extra comma makes sure the passed constant is a tuple as required
# Once variable created, the file will reject repeated creation
try:
    varmmd = grp.createVariable("mmd","f8",('timeutc',))
except RuntimeError:
    varmmd = grp['mmd']
varmmd[:] = tmpmmd
file.close()


#%%

# Destructively concatenate RASTA data to the merged dataset
outputfile = 'tmp/darwin_rasta.nc'
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

rastafilenames = glob.glob('/data/mcfarq/a/szhu28/research/HIWC/data/fulldataDarwin/RASTA/data/*.nc')

rastaset = [ read_rasta(filepath) for filepath in rastafilenames ]
timeutc = xr.open_dataset('tmp/darwin/darwin_sync_bulk.nc',group='/saffire').timeutc

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
dirpaths = glob.glob('/data/gpm/a/shared/szhu28/hiwcproc/darwin/mtsat/*')
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

sf = xr.open_dataset('tmp/darwin/darwin_sync_bulk.nc',group='/saffire')[ ['latitude','longitude'] ]
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
 'reflectance_nir',
 'temperature_sir',
 'temperature_67',
 'temperature_ir',
 'temperature_sw',
 'broadband_shortwave_albedo',
 'broadband_longwave_flux',
 'ir_cloud_emittance',
 'cloud_phase',
 'visible_optical_depth',
 'particle_size',
 'liquid_water_path',
 'cloud_effective_temperature',
 'cloud_top_pressure',
 'cloud_effective_pressure',
 'cloud_bottom_pressure',
 'cloud_top_height',
 'cloud_effective_height',
 'cloud_bottom_height']

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
outputfile = 'tmp/darwin_mtsat.nc'
warnings.filterwarnings('ignore')
for i in range(len(indarr)-1):
    print(i)
# for i in range(2):
    indi = slice(indarr[i],indarr[i+1])
    filepath = str(satfiles.sel(timesat=inserted[indarr[i]]).values)

    # sadly in the original data file, the missing_value attribute is a string instead of a number, 
    # and thus the auto NAN identification of xarray doesn't work
    tmpds = xr.open_dataset(filepath, mask_and_scale=False)
    tmpds = tmpds[varkeys]

    # longitude and latitude have different missing values
    for key in varkeys[0:2]:
        tmpds[key] = tmpds[key].where(~(tmpds[key]==-99999.))
    for key in varkeys[2:]:
        tmpds[key] = tmpds[key].where(~(tmpds[key]==-9.))

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
tmpmtsatproc['validlagmask'] = xr.DataArray( tmp <= np.timedelta64(30,'m') , dims=['timeutc'], coords={'timeutc':timeutc})
mtsatproc = xr.merge([tmpmtsatproc,mtsatproc],join='left')
mtsatproc.to_netcdf(outputfile,mode='w',
    format='NETCDF4',engine='netcdf4',group='/mtsatproc',unlimited_dims=['timeutc'])


#%%

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

