"""
Satellite Verification package
v 1.0

Framework for Class-1 verification (non-ice surface fields)

SST GHRSST (NCEI from PO.DAAC 1-day)
SSS SMAP (from PO.DAAC (8-day centered mean))
SSS SMOS (from NESDIS (3-day right side mean))
SSH MADT (from Copernicus (6-day left side mean))

Uses WOA13 climatology (SST, SSS) and 30-yr SSH avg from HYCOM reanalysis
Gaussian-weighted nearest neighbor interpolation for regridding, linear interp 
for climo and ice interpolations.

22 subregions (usual Global RTOFS subregions)
All forecasts (daily 0-192 hrs)

Notes:
    v0.1: this is a data-centric set of routines, where the validDate is
          derived from the observations, not the model.  All of the data 
          read routines need to support vDate=None to get the latest 
          data set.  Then the code will load the necessary model forecasts
          to match the data valid time.    
    v0.2: this is the MP version
    v0.3: split dbase into parameter tables to prevent collisions in r/w
    v0.4: cron-ready, added timeouts to sqlite3 to prevent db locking
    v0.5: switched smap from thredds to opendap server, added time attribute
          to all obs data sets (obs.attrs['time']) to standardize it.
    v0.6: fixed ssl cert issue with the po.daac opendap server, so smap data
          for the latest day can be determined.
    v0.7: switched from NCEP ice fields to CICE fields (from RTOFS) for ice
          masking.
    v0.8: ported to Mars
    v0.9: switched SMAP from RSS to JPL
    v1.0: reset SMAP, SMOS urls at NESDIS, GHRSST to NCEI at PO.DAAC, new interp
          for climo, ice.  Using only one processor for NN interp (regrid)
    v1.1: modified for RTOFSv2 nc files

    Todd Spindler

    v1.2: ported to Chen's Hera account with minor fixes, 1 March 2022    
"""
import matplotlib
matplotlib.use('agg')
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib.dates import MonthLocator, DateFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from mmab_toolkit import Basemap
import pyresample as pyr
import sqlite3
import pandas as pd
from datetime import datetime, timedelta
import concurrent.futures as cf
from mmab_toolkit import MidpointNormalize
from mmab_toolkit import add_mmab_logos
from skimage.measure import compare_ssim
from sklearn.metrics import mean_squared_error
import io
from PIL import Image
from glob import glob
import warnings
import os, sys
#import ipdb   # do not use this when running on compute farm!

#xr.set_options(keep_attrs=True)

WANT_POOL=True
maxjobs=12
USE_OSTIA_ICE=True
PLOT_FIELDS=True
UPDATE_DB=True
PLOT_STATS=True
MOVE_TO_POLAR=False
CLEANUP_IMAGES=False

srcDir='/scratch2/NCEPDEV/ovp/Lichuan.Chen/VPPPG/Global_RTOFS/EMC_ocean-verification/satellite'
refDir=srcDir+'/fix'
dbfile=refDir+'/satellite.db'

warnings.filterwarnings("ignore")

#-------------------------
def init_db(dbfile,table):
    """
    initialize the database, pass in the table name
    """
    # connect to db file
    conn = sqlite3.connect(dbfile,detect_types=sqlite3.PARSE_DECLTYPES,timeout=30.0)
    c = conn.cursor()    
    
    # Create table
    c.execute('CREATE TABLE IF NOT EXISTS '+table+' '+
              '(vdate timestamp,' +
              'fcst int,' +
              'platform string,' +
              'param string,' +
              'region string,' +
              'points int,' +
              'mean real,' +
              'std real,' + 
              'anom_mean real,' +
              'anom_std real,' + 
              'diff_mean real,' +
              'diff_std real,' + 
              'rmse real,' +               
              'bias real,' + 
              'cc real,' + 
              'si real,' + 
              'ssim real,' + 
              'anom_ssim real,' + 
              'unique(vdate,fcst,param,region))')
                                  
    # Save (commit) the changes
    conn.commit()
                 
    # close the connection
    conn.close()    
    return
#-------------------------
def update_db(dbfile,vdate,fcst,platform,param,region,points,mean,std,
              anom_mean,anom_std,diff_mean,diff_std,rmse,bias,cc,si,ssim,anom_ssim):
    """
    manage the SQLite database
    """
    
    init_db(dbfile,platform)  # if needed
        
    # connect to db file
    conn = sqlite3.connect(dbfile,detect_types=sqlite3.PARSE_DECLTYPES,timeout=30.0)
    c = conn.cursor()    
    
    # Larger example that inserts many records at a time
    update=(vdate,fcst,platform,param,region,points,mean,std,
            anom_mean,anom_std,diff_mean,diff_std,rmse,bias,cc,si,ssim,anom_ssim)

    c.execute(('REPLACE INTO '+platform+' VALUES ('+'?,'*len(update))[:-1]+')', update)
    
    # Save (commit) the changes
    conn.commit()
                 
    # close the connection
    conn.close()    
    
    return
#---------------------------
def read_db(dbfile,platform):
    """
    read the database back in 
    """
    conn=sqlite3.connect(dbfile,detect_types=sqlite3.PARSE_DECLTYPES,timeout=30.0)
    results=pd.read_sql('SELECT DISTINCT * FROM '+platform,conn)
    conn.close()
    
    return results
#-------------------------------------
# figure optimization and compression 
#-------------------------------------
def saveImage(imagefile,**kwargs):
    ram = io.BytesIO()
    plt.gcf().savefig(ram,**kwargs)
    ram.seek(0)
    im=Image.open(ram)
    im2 = im.convert('RGB').convert('P', palette=Image.WEB)
    im.save(imagefile, format='PNG',optimize=True)
    ram.close()
    return
#------------------------------------------------------------------------
def initialize_regions():
    """
    Standard Global RTOFS subregions
    """
    regions = {}  # python 3.6+ dicts preserve order 
    regions['global'] = {'name':'global','longname':'Global','lat':(-80,90),'lon':(0,360)}
    regions['natl'] = {'name':'natl','longname':'North Atlantic','lat':(0,60),'lon':(260,390)}
    regions['satl'] = {'name':'satl','longname':'South Atlantic','lat':(-80,0),'lon':(280,390)}
    regions['eqatl'] = {'name':'eqatl','longname':'Equatorial Atlantic','lat':(-30,30),'lon':(280,390)}
    regions['npac'] = {'name':'npac','longname':'North Pacific','lat':(0,70),'lon':(95,290)}
    regions['spac'] = {'name':'spac','longname':'South Pacific','lat':(-80,0),'lon':(95,295)}
    regions['eqpac'] = {'name':'eqpac','longname':'Equatorial Pacific','lat':(-30,30),'lon':(95,295)}
    regions['southern'] = {'name':'southern','longname':'Southern Ocean','lat':(-80,-30),'lon':(0,360)}
    regions['indian'] = {'name':'indian','longname':'Indian Ocean','lat':(-75,30),'lon':(10,130)}
    regions['arctic'] = {'name':'arctic','longname':'Arctic Ocean','lat':(47,90),'lon':(0,360)}
    regions['med'] = {'name':'med','longname':'Mediterranean Sea','lat':(29,48),'lon':(350,405)}
    regions['gulfstream'] = {'name':'gulfstream','longname':'Gulf Stream Region','lat':(15,47),'lon':(260,306)}
    regions['kuroshio'] = {'name':'kuroshio','longname':'Kuroshio Region','lat':(15,45),'lon':(115,155)}
    regions['agulhas'] = {'name':'agulhas','longname':'Agulhas Region','lat':(-50,-15),'lon':(0,45)}
    regions['somalia'] = {'name':'somalia','longname':'Somalia Region','lat':(-35,15),'lon':(25,60)}
    regions['northbrazil'] = {'name':'northbrazil','longname':'North Brazil Region','lat':(-15,15),'lon':(295,340)}
    regions['persiangulf'] = {'name':'persiangulf','longname':'Somalia & Persian Gulf','lat':(0,35),'lon':(40,65)}
    regions['australia'] = {'name':'australia','longname':'Australia & New Zealand','lat':(-50,-10),'lon':(110,180)}
    regions['indonesia'] = {'name':'indonesia','longname':'Indonesia & Philippines','lat':(-10,25),'lon':(100,130)}
    regions['gulfmex'] = {'name':'gulfmex','longname':'Gulf of Mexico','lat':(15,35),'lon':(260,285)}
    regions['westconus'] = {'name':'westconus','longname':'Western CONUS','lat':(10,71),'lon':(202,282)}
    regions['alaska'] = {'name':'alaska','longname':'Alaska','lat':(46,75),'lon':(180,218)}
    
    return regions
#---------------------------
def get_model(vDate,platform,param):
    """
    read RTOFS data and aggregate as needed by obs platforms
    aviso ssh = 1-day composite
    smap sss = 1-day composite
    smos sss = 1-day composite
    
    Load the correct files and average them in time.    
    Also note that obs are all at 12Z and model is at 0Z
    
    """
    catalog={'sss':'prog','sst':'prog','ssh':'diag','ice_coverage':'ice'}
    units={'sss':'PSU','sst':'deg C','ssh':'meters','ice_coverage':'fraction'}
         
    dataDir='/scratch2/NCEPDEV/ovp/Lichuan.Chen/Global/archive'
    refDir='/scratch2/NCEPDEV/ovp/Lichuan.Chen/VPPPG/Global_RTOFS/EMC_ocean-verification/satellite/fix'
    
    print(' ')
    for fcst in np.arange(0,193,24):
        agg={'smap':pd.date_range(vDate,vDate)-timedelta(fcst/24.),
             'smos':pd.date_range(vDate,vDate)-timedelta(fcst/24.),
             'aviso':pd.date_range(vDate,vDate)-timedelta(fcst/24.),
             'ghrsst':pd.date_range(vDate,vDate)-timedelta(fcst/24.)}
             
        print('processing',vDate,fcst)
        
        # here's where things get complicated
        # create a list of nc files with correct time span and assorted other bits
                
        template='{}/{}/rtofs_glo_2ds_{}{:03n}_{}.nc'
        
        if fcst > 0:
            ftype='f'
            fcst_hrs=fcst
        else:
            ftype='n'
            fcst_hrs=24
        
        filelist=[template.format(dataDir,i.strftime('%Y%m%d'),ftype,fcst_hrs,catalog[param]) for i in agg[platform]]
                    
        indata=xr.open_mfdataset(filelist,decode_times=True)
        
        indata=indata.mean(dim='MT')  # aggregation step
        #indata[param]=indata[param].astype('single')
        indata = indata[param][:-1,]  # select the parameter and fix the bad lons issue
        indata.coords['time']=vDate
        indata.coords['fcst']=fcst
        #indata.attrs['filelist']={fcst:filelist}
        try:
            outdata=xr.concat((outdata,indata),dim='fcst',coords='all')
        except:  # first go round
            outdata=indata.copy()            
            
    outdata=outdata.rename({'Longitude':'lon','Latitude':'lat',})
    # all coords need to be single precision
    outdata['lon']=outdata.lon.astype('single')
    outdata['lat']=outdata.lat.astype('single')
    outdata.attrs['platform']='rtofs '+platform
    outdata.attrs['units']=units[param]
    return outdata

#------------------------------------------------------------------
def get_ice(vDate=None):

    if vDate==None:
        files=glob('/scratch2/NCEPDEV/ovp/Lichuan.Chen/OSTIA/OSTIA*.nc')
        files.sort()
        icefile=files[-1]
    else:
        icefile='/scratch2/NCEPDEV/ovp/Lichuan.Chen/OSTIA/OSTIA-UKMO-L4-GLOB-v2.0_'+vDate.strftime('%Y%m%d')+'.nc'
        
    # load ice data
    print('reading OSTIA ice')
    
    if not os.path.exists(icefile):
        print('missing OSTIA ice file for',vDate)
        return
    
    data=xr.open_dataset(icefile,decode_times=True)
    data=data.rename({'sea_ice_fraction':'ice'})
    
    # all coords need to be single precision
    data2=data.ice.astype('single')
    data2['lon']=data2.lon.astype('single')
    data2['lat']=data2.lat.astype('single')
        
    return data2
#-------------------------------------------------------
def get_ghrsst(vDate=None):
    """
    NESDIS GHRSST OSPO product daily day/night
    daily composite
    """    
    if vDate==None:
        files=glob('/scratch2/NCEPDEV/ovp/Lichuan.Chen/GHRSST/GHRSST*.nc')
        files.sort()
        url=files[-1]    
    else:
        url='/scratch2/NCEPDEV/ovp/Lichuan.Chen/GHRSST/GHRSST-OSPO-L4-GLOB_'+vDate.strftime('%Y%m%d')+'.nc'
        
    data=xr.open_dataset(url,decode_times=True)
    data['time']=data.time-pd.Timedelta('12H')  # shift 12Z offset time to 00Z
    data2=data.analysed_sst.astype('single')-273.15 # convert from Kelvin
            
    print('Retrieved GHRSST data from NESDIS for',data2.time.values)
            
    #data2=data2.rename({'longitude':'lon','latitude':'lat'})        
    
    # all coords need to be single precision
    data2['lon']=data2.lon.astype('single')
    data2['lat']=data2.lat.astype('single')
    data2.attrs['platform']='ghrsst'
    data2.attrs['units']='Celsius'
    data2.attrs['time']=pd.Timestamp(data2.time.values[0].strftime())
    
    return data2
#-------------------------------------------------------
def get_smos(vDate=None):
    """
    NESDIS SMOS product daily aggregate Level 3 grid      
    daily composite
    """
    if vDate==None:
        files=glob('/scratch2/NCEPDEV/ovp/Lichuan.Chen/SMOS/SMOS*.nc')
        files.sort()
        url=files[-1]    
    else:
        url='/scratch2/NCEPDEV/ovp/Lichuan.Chen/SMOS/SMOS-L3-GLOB_'+vDate.strftime('%Y%m%d')+'.nc'

    data=xr.open_dataset(url,decode_times=True)
    data['time']=data.time-pd.Timedelta('12H')  # shift 12Z offset time to 00Z
    data2=data['sss'].astype('single')
            
    print('Retrieved SMOS data from NESDIS for',data2.time.values)
            
    data2=data2.rename({'longitude':'lon','latitude':'lat'})        
    
    # all coords need to be single precision
    data2['lon']=data2.lon.astype('single')
    data2['lat']=data2.lat.astype('single')
    data2.attrs['platform']='smos'
    data2.attrs['units']='PSU'
    data2.attrs['time']=pd.Timestamp(data2.time.values[0])
    
    return data2
#-------------------------------------------------------
def get_smap(vDate=None):
    """
    NESDIS SMAP product daily aggregate Level 3 grid      
    daily composite
    """    
    if vDate==None:
        files=glob('/scratch2/NCEPDEV/ovp/Lichuan.Chen/SMAP/SMAP*.nc')
        files.sort()
        url=files[-1]    
    else:
        url='/scratch2/NCEPDEV/ovp/Lichuan.Chen/SMAP/SMAP-L3-GLOB_'+vDate.strftime('%Y%m%d')+'.nc'
    
    data=xr.open_dataset(url,decode_times=True)
    data['time']=data.time-pd.Timedelta('12H')  # shift 12Z offset to 00Z    
    data2=data['sss'].astype('single')
                    
    print('Retrieved SMAP data from NESDIS for',data2.time.values)
                
    # all coords need to be single precision
    data2['lon']=data2.lon.astype('single')
    data2['lat']=data2.lat.astype('single')
    data2.attrs['platform']='smap'
    data2.attrs['units']='PSU'
    data2.attrs['time']=pd.Timestamp(data2.time.values[0]) 
    
    return data2
#---------------------------------------------------------------
def get_aviso(vDate=None):
    """
    AVISO is hosted on Copernicus CMEMS
    ssh/sla data is double precision on the server, so convert to single
    """
    
    if vDate==None:
        files=glob('/scratch2/NCEPDEV/ovp/Lichuan.Chen/AVISO/nrt_*.nc')
        files.sort()
        filename=files[-1]    
    else:
        filename='/scratch2/NCEPDEV/ovp/Lichuan.Chen/AVISO/nrt_global_allsat_phy_l4_'+vDate.strftime('%Y%m%d')+'.nc'
                
    data=xr.open_dataset(filename,decode_times=True)
    
    sla=data.sla.astype('single')
    sla.attrs['platform']='aviso'
    sla.attrs['filename']=filename.split('/')[-1]
    sla.attrs['time']=pd.Timestamp(data.time.values[0])
    sla=sla.rename({'longitude':'lon','latitude':'lat'})
    # all coords need to be single precision
    sla['lon']=sla.lon.astype('single')
    sla['lat']=sla.lat.astype('single')
    sla.attrs['units']='meters'

    adt=data.adt.astype('single')
    adt.attrs['platform']='aviso'
    adt.attrs['filename']=filename
    adt.attrs['time']=pd.Timestamp(data.time.values[0])
    adt=adt.rename({'longitude':'lon','latitude':'lat'})
    # all coords need to be single precision
    adt['lon']=adt.lon.astype('single')
    adt['lat']=adt.lat.astype('single')
    adt.attrs['units']='meters'
    
    return sla, adt
#------------------------------------------------
def get_woa(vDate,param):
    """
    Read World Ocean Atlas climatology data        
    WOA monthly data is centered on mid-month, need to scale the vDate against
    the date range to create a weighted average of the straddling month fields.
    """    
    climoDir='/scratch2/NCEPDEV/ovp/Lichuan.Chen/Global/climo/WOA'
    vDate=pd.Timestamp(vDate) # just in case
    if vDate.day==15:  # even for Feb, just because
        if param=='sst':
            climofile="woa13_decav_t{:02n}_04v2.nc".format(vDate.month)    
            data=xr.open_dataset(climoDir+'/'+climofile,decode_times=False)
            data=data['t_an'].squeeze()[0,]  # surface only
        elif param=='sss':
            climofile="woa13_decav_s{:02n}_04v2.nc".format(vDate.month)    
            data=xr.open_dataset(climoDir+'/'+climofile,decode_times=False)
            data=data['s_an'].squeeze()[0,]  # surface only
    else:  # need to scale things
        if vDate.day < 15:
            start=pd.Timestamp(vDate.year,vDate.month,15)+pd.tseries.offsets.DateOffset(months=-1)
            stop=pd.Timestamp(vDate.year,vDate.month,15)
        else:
            start=pd.Timestamp(vDate.year,vDate.month,15)
            stop=pd.Timestamp(vDate.year,vDate.month,15)+pd.tseries.offsets.DateOffset(months=1)            
        left=(vDate-start)/(stop-start)
        #right=(stop-vDate)/(stop-start)
        if param=='sst':
            climofile1="woa13_decav_t{:02n}_04v2.nc".format(start.month)    
            climofile2="woa13_decav_t{:02n}_04v2.nc".format(stop.month)    
            data1=xr.open_dataset(climoDir+'/'+climofile1,decode_times=False)
            data2=xr.open_dataset(climoDir+'/'+climofile2,decode_times=False)
            data1=data1['t_an'].squeeze()[0,]  # surface only
            data2=data2['t_an'].squeeze()[0,]  # surface only
        elif param=='sss':
            climofile1="woa13_decav_s{:02n}_04v2.nc".format(start.month)    
            climofile2="woa13_decav_s{:02n}_04v2.nc".format(stop.month)    
            data1=xr.open_dataset(climoDir+'/'+climofile1,decode_times=False)
            data2=xr.open_dataset(climoDir+'/'+climofile2,decode_times=False)
            data1=data1['s_an'].squeeze()[0,]  # surface only
            data2=data2['s_an'].squeeze()[0,]  # surface only
            
        #data=data1*left+data2*right
        data=data1+((data2-data1)*left)
        climofile='weighted average of '+climofile1+' and '+climofile2
    # all coords need to be single precision
    data['lon']=data.lon.astype('single')
    data['lat']=data.lat.astype('single')
    data.attrs['platform']='woa'
    data.attrs['filename']=climofile

    return data
#------------------------------------------------
def get_hycom_climo(vDate):
    """
    Read HYCOM 1995-2015 mean data
    This is used for the mssh values
    HYCOM monthly data is centered on mid-month, need to scale the vDate against
    the date range to create a weighted average of the straddling month fields.
    """
    climoDir='/scratch2/NCEPDEV/ovp/Lichuan.Chen/Global/climo/HYCOM'
    vDate=pd.Timestamp(vDate) # just in case
    if vDate.day==15:  # even for Feb, just because
        climofile="hycom_GLBv0.08_53X_archMN.1994_{0:02n}_2015_{0:02n}_ssh.nc".format(vDate.month)
        data=xr.open_dataset(climoDir+'/'+climofile,decode_times=False)
        data=data['surf_el'].copy().squeeze()
    else:  # need to scale things
        if vDate.day < 15:
            start=pd.Timestamp(vDate.year,vDate.month,15)+pd.tseries.offsets.DateOffset(months=-1)
            stop=pd.Timestamp(vDate.year,vDate.month,15)
        else:
            start=pd.Timestamp(vDate.year,vDate.month,15)
            stop=pd.Timestamp(vDate.year,vDate.month,15)+pd.tseries.offsets.DateOffset(months=1)
        left=(vDate-start)/(stop-start)
        #right=(stop-vDate)/(stop-start)
        climofile1="hycom_GLBv0.08_53X_archMN.1994_{0:02n}_2015_{0:02n}_ssh.nc".format(start.month)
        climofile2="hycom_GLBv0.08_53X_archMN.1994_{0:02n}_2015_{0:02n}_ssh.nc".format(stop.month)
        data1=xr.open_dataset(climoDir+'/'+climofile1,decode_times=False)
        data2=xr.open_dataset(climoDir+'/'+climofile2,decode_times=False)
        data1=data1['surf_el'].copy().squeeze()
        data2=data2['surf_el'].copy().squeeze()                    
        #data=data1*left+data2*right
        data=data1+((data2-data1)*left)
        climofile='weighted average of '+climofile1+' and '+climofile2        
    data.coords['time']=datetime(vDate.year,vDate.month,1)   # just a reference to the month 
    # all coords need to be single precision
    data['lon']=data.lon.astype('single')
    data['lat']=data.lat.astype('single')
    data.attrs['platform']='hycom'
    data.attrs['filename']=climofile
    
    return data
#-------------------------------------------------
def regrid(model,obs):
    """
    regrid data to obs -- this assumes DataArrays   
    """
    model2=model.copy()
    model2_lon=model2.lon.values
    model2_lat=model2.lat.values
    model2_data=model2.to_masked_array()        
    if model2_lon.ndim==1:
        model2_lon,model2_lat=np.meshgrid(model2_lon,model2_lat)
       
    obs2=obs.copy()
    obs2_lon=obs2.lon.astype('single').values
    obs2_lat=obs2.lat.astype('single').values
    obs2_data=obs2.astype('single').to_masked_array()    
    if obs2.lon.ndim==1:
        obs2_lon,obs2_lat=np.meshgrid(obs2.lon.values,obs2.lat.values)
        
    model2_lon1=pyr.utils.wrap_longitudes(model2_lon)
    model2_lat1=model2_lat.copy()
    obs2_lon1=pyr.utils.wrap_longitudes(obs2_lon)
    obs2_lat1=obs2_lat.copy()
    
    # pyresample gausssian-weighted kd-tree interp 
    # define the grids
    orig_def = pyr.geometry.GridDefinition(lons=model2_lon1,lats=model2_lat1)
    targ_def = pyr.geometry.GridDefinition(lons=obs2_lon1,lats=obs2_lat1)
    radius=50000
    sigmas=25000
    model2_data2=pyr.kd_tree.resample_gauss(orig_def,model2_data,targ_def,
                                            radius_of_influence=radius,
                                            sigmas=sigmas,
                                            fill_value=None)        
    # match the masks
    #combined_mask=np.ma.mask_or(model2_data2.mask,obs2_data.mask)
    #model2_data2.mask=combined_mask
    #obs2_data.mask=combined_mask
        
    model=xr.DataArray(model2_data2,coords=[obs.lat.values,obs.lon.values],dims=['lat','lon'])
    #obs=xr.DataArray(obs2_data,coords=[obs2.lat.values,obs2.lon.values],dims=['lat','lon'])
    
    return model

#---------------------------
def calc_stats(obs,model,climo,param,obs_anom=None):
    """
    compute anomaly stats     
    """
    stats={}
    model_anom=model-climo
    if param != 'ssh':
        obs_anom=obs-climo
    diff=model-obs
    bias = diff.mean()
    obs_sdev=obs_anom.std()
    model_sdev=model_anom.std()
    if (obs.compressed().size > 0) and (model.compressed().size > 0): 
        rmse=mean_squared_error(obs.compressed(),model.compressed())**0.5
    else:
        rmse=np.nan
    cc=np.ma.corrcoef(obs_anom.compressed(),model_anom.compressed())[0,1]
    si=100.0*((diff**2).mean()**0.5 - bias**2)/obs.mean()
    
    # self similarity index metric computed on non-masked areas of the model and obs anomaly grids
    #data_range=max((obs.max(),model.max()))-min((obs.min(),model.min()))
    #anom_data_range=max((obs_anom.max(),model_anom.max()))-min((obs_anom.min(),model_anom.min()))
    ssim = compare_ssim(obs.filled(fill_value=999).astype('single'), 
                        model.filled(fill_value=999).astype('single'))
    ssim_anom = compare_ssim(obs_anom.filled(fill_value=999).astype('single'), 
                             model_anom.filled(fill_value=999).astype('single'))
    
    stats['mean']=model.mean()
    stats['std']=model.std()
    stats['anom_mean']=model_anom.mean()
    stats['anom_std']=model_anom.std()    
    stats['diff_mean']=diff.mean()
    stats['diff_std']=diff.std()
    stats['rmse']=rmse
    stats['bias']=bias
    stats['cc']=cc
    stats['si']=si
    stats['ssim']=ssim
    stats['anom_ssim']=ssim_anom
    stats['points']=obs.compressed().size
    
    return stats
#-----------------------------------------------
def plot_stats(vDate,platform,param,region,imageDir=None,lookback='60 d'):
    
    # read in the database
    df=read_db(dbfile,platform)
    df=df[df.region==region['name']]
    
    df.index=df.vdate
    minDate=pd.to_datetime(vDate)-pd.Timedelta(lookback)
    df=df[minDate:vDate]
    
    #df.index=df.vdate     
    #minDate=pd.to_datetime(vDate)-pd.Timedelta(lookback)
    #df=df[minDate:vDate]

    # stats plots
    stats={'anom_mean':'Anomaly Mean',
           'rmse':'RMSE',
           'bias':'Bias',
           'cc':'Correlation Coeff',
           'si':'Scatter Index',
           'anom_ssim':'Anomaly Self-Similarity Index'}
    
    for stat,statname in stats.items():
        for fcst in range(24,193,24):
            fig=plt.figure(dpi=150)
            ax=fig.add_subplot(111)
            df2=df[(df.fcst==fcst)].copy()
            df2=df2.asfreq('D') # create continuous timeseries with gaps
            matplotlib.rcParams['legend.fontsize'] = 'x-small'
            matplotlib.rcParams['legend.facecolor'] = 'white'
            # label all fcst lines (green) using 24 hr fcst for the single label
            for fcst2 in range(24,193,24):
                if fcst2==24:
                    #df[df.fcst==fcst2].plot(ax=ax,x='vdate',y=stat,x_compat=True,color='green',alpha=0.5,label='all forecasts',fontsize='x-small')
                    df[df.fcst==fcst2].asfreq('D').plot(ax=ax,x='vdate',y=stat,x_compat=True,
                    	color='green',linewidth=0.5,label='all forecasts',fontsize='x-small')
                else:
                    #df[df.fcst==fcst2].plot(ax=ax,x='vdate',y=stat,x_compat=True,color='green',alpha=0.5,label='_none',fontsize='x-small')
                    ##df[df.fcst==fcst2].plot(ax=ax,x='vdate',y=stat,x_compat=True,color='green',linewidth=0.5,label='_none',fontsize='x-small')
                    df[df.fcst==fcst2].asfreq('D').plot(ax=ax,x='vdate',y=stat,x_compat=True,
                    	color='green',linewidth=0.5,fontsize='x-small',label='_none')
            # label nowcast line (blue dashed)
            df[df.fcst==0].asfreq('D').plot(ax=ax,x='vdate',grid=True,y=stat,x_compat=True,
            	color='blue',linestyle='--',label='nowcast',fontsize='x-small')
            # label fcst line (red solid)
            df2.plot(ax=ax,x='vdate',y=stat,grid=True,x_compat=True,
            	color='red',label="{:03n} hr fcst".format(fcst),fontsize='x-small')
            # add number of points
            ax2=df2.plot(ax=ax,x='vdate',y='points',grid=False,secondary_y=True,sharex=True,x_compat=True,
            	linestyle='--',color='gray',label='number of points',fontsize='x-small')
            ax2.ticklabel_format(axis='y',style='sci',useOffset=False)
            if np.sign(np.prod(ax.get_ylim())) <= 0:
                ax.axhline(color='k',linestyle='-',linewidth=1)

            # set up the single legend for both axes
            lines_1, labels_1 = ax.get_legend_handles_labels()
            lines_2, labels_2 = ax2.get_legend_handles_labels()

            lines = lines_1 + lines_2
            labels = labels_1 + labels_2
            ax.legend(lines, labels, framealpha=0.8)
            ax.grid()
            ax.set_xlabel('Valid Date',fontsize='small')
            ax.set_ylabel(statname,fontsize='small')
            ax2.set_ylabel('Number of points',fontsize='small')
            title=(param.upper()+' '+statname+' '+region['longname']+' '+vDate+' '+str(fcst)).replace('SSH','SLA')
            plt.title(title,fontsize='medium')
            add_mmab_logos()
            if not os.path.isdir(imageDir):
                os.makedirs(imageDir)
            #imageFile=('class1_stats_'+param+'_'+stat+'_'+region['name']+'_'+vDate+'_'+str(fcst)+'.png').replace('ssh','sla')
            imageFile="class1_stats_{}_{}_{}_{}_{}_{:03n}.png".format(platform,param.replace(' ','_'),
                                                                   stat,
                                                                   region['name'],
                                                                   vDate,
                                                                   fcst)
            plt.savefig(imageDir+'/'+imageFile)
            #saveImage(imageDir+'/'+imageFile)
            plt.close()
        
    return
#-----------------
# cartopy version 
#-----------------
def make_map_cartopy(data,platform,region,vDate,fcst,param,imageDir=None,cmap=None):
    """
    cartopy graphics
    """
    print('mapping',param,platform.lower(),region['name'],vDate,fcst)
    # set colormap as needed
    if np.sign(data.min())*np.sign(data.max())==-1 and cmap==None:
        cmap='bwr'    
    else:
        cmap='gist_rainbow_r'
    #data['lon']=np.remainder(data.lon+180,360)-180  # shift lat/lon to -180/180
    plt.close()
    plt.figure(dpi=150)
    if region['name']=='arctic':
        crs=ccrs.NorthPolarStereo()
    elif region['name']=='southern':
        crs=ccrs.SouthPolarStereo()
    else:
        #crs=ccrs.Mercator(central_longitude=data.lon.values.mean())
        crs=ccrs.PlateCarree(central_longitude=data.lon.values.mean())
        
    #ax=plt.axes(projection=ccrs.Mercator(central_longitude=np.remainder(data.lon.mean()+180,360)-180))
    ax=plt.subplot(111,projection=crs)

    if region['name'] != 'global' and region['name'] != 'arctic' and region['name'] != 'southern':
        extent=region['lon']+region['lat']
        #ax.set_extent(extent,crs=ccrs.PlateCarree())
        ax.set_extent(extent)
    else:
        # Compute a circle in axes coordinates, which we can use as a boundary
        # for the map. We can pan/zoom as much as we like - the boundary will be
        # permanently circular.
        if region['name'] == 'arctic' or region['name'] == 'southern':
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)        
            ax.set_boundary(circle, transform=ax.transAxes)
        
    ax.add_feature(cfeature.LAND,zorder=1)  # eventually
    if region['name'] != 'global':
        ax.coastlines(resolution='10m',zorder=1)
    else:
        ax.coastlines(zorder=1)
    if region['name'] == 'arctic' or region['name'] == 'southern':
        gl=ax.gridlines()
    else:
        gl=ax.gridlines(draw_labels=True)
        gl.ylabels_right=False
        gl.xlabels_top=False

    data.plot(x='lon',y='lat',transform=ccrs.PlateCarree(),
    	vmin=data.min(),vmax=data.max(),
        cbar_kwargs={'shrink':0.75},
        #cbar_kwargs={'orientation':'horizontal','ax':ax,'use_gridspec':True},
        cmap=cmap,
        norm=MidpointNormalize(midpoint=0.,vmin=data.min(),vmax=data.max()))
    
    ax.xaxis.set_tick_params(labelsize='x-small')
    ax.yaxis.set_tick_params(labelsize='x-small')
                   
    title=param.upper()+' '+platform.upper()+' '+region['longname']+' '+vDate+' fcst: '+str(fcst)
    plt.title(title.replace('SMOS SMOS','SMOS SSS').replace('SMAP SMAP','SMAP SSS').replace('GHRSST GHRSST','GHRSST SST'),fontsize='small')    
    
    add_mmab_logos()
    imagename="class1_{}_{}_{}_{}_{:03n}.png".format(platform,
                                                  param.replace(' ','_').lower(),
                                                  vDate,
                                                  region['name'],
                                                  fcst)
    #imagename=imagename.replace('smos_smos','smos_sss').replace('smap_smap','smap_sss')
    if not os.path.isdir(imageDir):
        os.makedirs(imageDir)
    plt.savefig(imageDir+'/'+imagename)
    #saveImage(imageDir+'/'+imagename)
    return
#-----------------
# basemap version 
#-----------------
def make_map_basemap(data,platform,region,vDate,fcst,param,units,limits=None,imageDir=None,cmap=None):
    """
    basemap graphics
    """
    
    print('mapping',platform,param,region['name'],vDate,fcst)
    # set colormap as needed
    if np.sign(data.min())*np.sign(data.max())==-1 and cmap==None:
        cmap='bwr'
    else:
        cmap='gist_rainbow_r'
    #data['lon']=np.remainder(data.lon+180,360)-180  # shift lat/lon to -180/180
    plt.close()
    plt.figure(dpi=150)
    
    if region['name']=='arctic':
        m=Basemap(projection='npstere',
                  boundinglat=data.lat.min(),lon_0=data.lon.mean(),
                  resolution='i',round=True)
    elif region['name']=='southern':
        m=Basemap(projection='spstere',
                  boundinglat=data.lat.max(),lon_0=data.lon.mean(),
                  resolution='i',round=True)
    else:
        m=Basemap(projection='mill',
                  llcrnrlon=data.lon.min(),urcrnrlon=data.lon.max(),
                  llcrnrlat=data.lat.min(),urcrnrlat=data.lat.max(),
                  resolution='i')
    #lon=data.lon.values
    #lat=data.lat.values
    #plot_data=data.to_masked_array()

    if limits==None:
        vmin=data.min()
        vmax=data.max()
        extend='neither'
    else:
        vmin=min(limits)
        vmax=max(limits)        
        if data.min() < vmin and data.max() > vmax:
            extend='both'
        else:
            if data.max() > vmax:
                extend='max'
            elif data.min() < vmin:
                extend='min'
            else:
                extend='neither'
            
    x,y=m(*np.meshgrid(data.lon.values,data.lat.values))

    if param.find('anomaly')!=-1 or param.find('sla')!=-1 or param.find('ssh')!=-1:
        m.pcolormesh(x,y,data,
                     vmin=vmin,vmax=vmax,
                     cmap=cmap,
                     norm=MidpointNormalize(midpoint=0.,
                     vmin=vmin,vmax=vmax))
    else:
        m.pcolormesh(x,y,data,
                     vmin=vmin,vmax=vmax,
                     cmap=cmap)

    plt.clim(vmin=vmin,vmax=vmax)
    cbar=m.colorbar(extend=extend,size='2.5%')
    cbar.ax.tick_params(labelsize='x-small')
    cbar.set_label(units)

    #m.drawlsmask(land_color='grey',ocean_color='none',lakes=True)
    m.fillcontinents()
    m.drawcoastlines(linewidth=0.2)
    m.gridlines()
    
    #title=param.upper()+' '+platform.upper()+' '+region['longname']+' '+vDate+' fcst: '+str(fcst)+ \
    #'\nMinimum: '+
    
    # only use forecast in model maps
    if platform.upper().find('RTOFS')!=-1:
        title="{} {} {} {} fcst {}\nMin: {:6.2f} Max: {:6.2f}".format(param.upper(),
                                                                      platform.upper(),
                                                                      region['longname'],
                                                                      vDate,
                                                                      fcst,
                                                                      data.min().values,
                                                                      data.max().values)
    else:
        title="{} {} {} {}\nMin: {:6.2f} Max: {:6.2f}".format(param.upper().replace('SSH ANOMALY','SLA'),
                                                              platform.upper(),
                                                              region['longname'],
                                                              vDate,
                                                              data.min().values,
                                                              data.max().values)
                                                                          
    plt.title(title,fontsize='small')
    
    add_mmab_logos()
    if platform.upper().find('RTOFS')!=-1:
        imagename="class1_{}_{}_{}_{}_{:03n}.png".format(platform.replace(' ','_'),
                                                         param.replace(' ','_').lower(),
                                                         vDate,
                                                         region['name'],
                                                         fcst)
    else:
        imagename="class1_{}_{}_{}_{}.png".format(platform,
                                                  param.replace(' ','_').lower(),
                                                  vDate,
                                                  region['name'])
        
    #imagename=imagename.replace('rtg_rtg','rtg_sst').replace('ghrsst_ghrsst','ghrsst_sst').replace('smos_smos','smos_sss').replace('smap_smap','smap_sss')
    if not os.path.isdir(imageDir):
        os.makedirs(imageDir)
    plt.savefig(imageDir+'/'+imagename,dpi=150)
    #saveImage(imageDir+'/'+imagename,dpi=150)
    return
#---------------------------
def expand_grid(data):
    """
    concatenate global data for edge wraps
    """
    
    data2=data.copy()
    data2['lon']=data2.lon+360
    data3=xr.concat((data,data2),dim='lon')
    return data3        
#-------------------------------
def process_data(model,obs,climo,ice,platform,param,fcst,imageDir,obs_anom=None):
    """
    process both model and obs
    use cice ice concentration as ice mask
    obs2_anom is sla from aviso
    """
    
    #ipdb.set_trace()
    
    # subregion can be computed right here
    regions=initialize_regions()

    print('regridding model to obs')
    model2=regrid(model,obs)

    # combine obs ice mask with ncep
    obs2=obs.to_masked_array()
    if param=='ssh':
        obs_anom2=obs_anom.to_masked_array()
    ice2=ice.to_masked_array()
    climo2=climo.to_masked_array()
    model2=model2.to_masked_array()
    
    # reconcile with obs 
    obs2.mask=np.ma.mask_or(obs2.mask,ice2>0.0)
    obs2.mask=np.ma.mask_or(obs2.mask,climo2.mask)
    obs2.mask=np.ma.mask_or(obs2.mask,model2.mask)
    climo2.mask=obs2.mask
    model2.mask=obs2.mask
    if param=='ssh':
        obs_anom2.mask=obs2.mask

    obs2=xr.DataArray(obs2,coords=[obs.lat.values,obs.lon.values], dims=['lat','lon'])
    if platform=='aviso':
        obs_anom2=xr.DataArray(obs_anom2,coords=[obs.lat.values,obs.lon.values], dims=['lat','lon'])
    model2=xr.DataArray(model2,coords=[obs.lat.values,obs.lon.values], dims=['lat','lon'])
    climo2=xr.DataArray(climo2,coords=[obs.lat.values,obs.lon.values], dims=['lat','lon'])
                    
    # now concat the data to handle edge cases (med, natl, satl)
    
    model2=expand_grid(model2)
    climo2=expand_grid(climo2)
    obs2=expand_grid(obs2)
    if platform=='aviso':
        obs_anom2=expand_grid(obs_anom2)

    if WANT_POOL:
        pool = cf.ProcessPoolExecutor(maxjobs)
        
    all_stats=[]
    for region in regions.values():    
        (latmin,latmax)=region['lat']
        (lonmin,lonmax)=region['lon']
        model3=model2.where((model2.lon>=lonmin)&(model2.lon<=lonmax)&
        	(model2.lat>=latmin)&(model2.lat<=latmax),drop=True)
        climo3=climo2.where((climo2.lon>=lonmin)&(climo2.lon<=lonmax)&
        	(climo2.lat>=latmin)&(climo2.lat<=latmax),drop=True)
        obs3=obs2.where((obs2.lon>=lonmin)&(obs2.lon<=lonmax)&
        	(obs2.lat>=latmin)&(obs2.lat<=latmax),drop=True)
        if platform=='aviso':
            obs_anom3=obs_anom2.where((obs_anom2.lon>=lonmin)&(obs_anom2.lon<=lonmax)&
            	(obs_anom2.lat>=latmin)&(obs_anom2.lat<=latmax),drop=True)
        #else:
        #    model3=model2
        #    climo3=climo2
        #    obs3=obs2
        #    if platform=='aviso':
        #        obs_anom3=obs_anom2
                               
        # compute statistics 
        if platform=='aviso':
            stats=calc_stats(obs3.to_masked_array(),model3.to_masked_array(),climo3.to_masked_array(),param,obs_anom3.to_masked_array())
        else:
            stats=calc_stats(obs3.to_masked_array(),model3.to_masked_array(),climo3.to_masked_array(),param)
        
        stats['vdate']=obs.attrs['time'].to_pydatetime()
        stats['fcst']=fcst
        stats['region']=region['name']
        all_stats.append(stats)
        
        # local valid date string for images 
        vDate=stats['vdate'].strftime('%Y%m%d')

        # map the data
        if PLOT_FIELDS:
            # choose basemap or cartopy here
            if platform=='aviso':
                cmap=None
            else:
                cmap='gist_rainbow_r'

            mdata=model3
            mcdata=model3-climo3
            odata=obs3
            if platform=='aviso':
                ocdata=obs_anom3
            else:
                ocdata=obs3-climo3
            
            # lock the colorbars
            mlim=[mdata.min(),mdata.max()]
            mclim=[mcdata.min(),mcdata.max()]
            olim=[odata.min(),odata.max()]
            oclim=[ocdata.min(),ocdata.max()]
            if param=='sss':
                mlim=[25,mdata.max()]
                olim=[25,odata.max()]
                
            if WANT_POOL:
                if fcst==24:
                    pool.submit(make_map_basemap,odata,obs.platform,region,vDate,fcst,param,obs.units,olim,imageDir,cmap)
                    pool.submit(make_map_basemap,ocdata,obs.platform,region,vDate,fcst,param+' anomaly',obs.units,oclim,imageDir)
                pool.submit(make_map_basemap,mdata,model.platform,region,vDate,fcst,param,obs.units,olim,imageDir,cmap)
                pool.submit(make_map_basemap,mcdata,model.platform,region,vDate,fcst,param+' anomaly',obs.units,oclim,imageDir)
            else:
                if fcst==24:
                    make_map_basemap(odata,obs.platform,region,vDate,fcst,param,obs.units,olim,imageDir,cmap)
                    make_map_basemap(ocdata,obs.platform,region,vDate,fcst,param+' anomaly',obs.units,oclim,imageDir)
                make_map_basemap(mdata,model.platform,region,vDate,fcst,param,obs.units,olim,imageDir,cmap)
                make_map_basemap(mcdata,model.platform,region,vDate,fcst,param+' anomaly',obs.units,oclim,imageDir)
                
    if WANT_POOL:
        print('shutting down pool')
        pool.shutdown(wait=True)
        
    return all_stats    
#---------------------------
# main routine              
#---------------------------
if __name__ == '__main__':
        
    platform=sys.argv[1]

    if len(sys.argv) == 3: 
        vDate=datetime.strptime(sys.argv[2],'%Y%m%d')
    else:
        vDate=None
        
    print('Starting Satellite V&V at',datetime.now(),'for',vDate)
        
    fcsts=dict(zip(range(9),range(0,193,24)))
    
    if platform=='aviso':
        param='ssh'
        obs_anom,obs=get_aviso(vDate)  # aviso return sla and adt
        obs_anom=obs_anom.squeeze()
        vDate=obs.attrs['time']
        model=get_model(vDate,platform,param)
        climo=get_hycom_climo(vDate)
    elif platform=='ghrsst':
        param='sst'
        obs=get_ghrsst(vDate)
        vDate=obs.attrs['time']
        model=get_model(vDate,platform,param) # RTG is a single time
        climo=get_woa(vDate,param)
    elif platform=='smap':
        param='sss'
        obs=get_smap(vDate)
        vDate=obs.attrs['time']
        model=get_model(vDate,platform,param)
        climo=get_woa(vDate,param)
    elif platform=='smos':
        param='sss'
        obs=get_smos(vDate)
        vDate=obs.attrs['time']
        model=get_model(vDate,platform,param)
        climo=get_woa(vDate,param)
    else:
        print('Platform',platform,'not recognized')
        sys.exit()
                
    obs=obs.squeeze()
                
    # need ice for masking.  Choose between CDAS and CICE-4 (model)
    if USE_OSTIA_ICE:
        ice=get_ice(vDate)
    else:
        ice=get_model(vDate,platform,'ice_coverage')
        ice=ice[0]  # only want nowcast ice

    # set the image directory here
    vDate=vDate.strftime('%Y%m%d')  # want it as a string 
    imageDir='/scratch2/NCEPDEV/stmp1/Lichuan.Chen/images/class-4/satellite/'+vDate

    print('regridding climo to obs')
    #climo2=regrid(climo,obs)
    climo=climo.squeeze()
    climo=regrid(climo,obs)
    print('\nregridding ice to obs')    
    #ice2=regrid(ice,obs)
    ice=regrid(ice,obs)

    for nfcst,fcst in fcsts.items():
        if platform == 'aviso':
            results=process_data(model[nfcst,],obs,climo,ice,platform,param,fcst,imageDir,obs_anom)
        else:
            results=process_data(model[nfcst,],obs,climo,ice,platform,param,fcst,imageDir)
        if UPDATE_DB:
            for result in results:
                update_db(dbfile,
                          result['vdate'],
                          fcst,
                          platform,
                          param,
                          result['region'],
                          result['points'],
                          result['mean'],
                          result['std'],
                          result['anom_mean'],
                          result['anom_std'],
                          result['diff_mean'],
                          result['diff_std'],
                          result['rmse'],
                          result['bias'],
                          result['cc'],
                          result['si'],
                          result['ssim'],
                          result['anom_ssim'])
            
    if PLOT_STATS:
        regions=initialize_regions()
                
        if WANT_POOL:
            pool = cf.ProcessPoolExecutor(maxjobs)
                
        for name,region in regions.items():
            if WANT_POOL:
                pool.submit(plot_stats,vDate,platform,param,region,imageDir)
            else:
                plot_stats(vDate,platform,param,region,imageDir)
                
        if WANT_POOL:
            print('shutting down pool')
            pool.shutdown(wait=True)
        
    if MOVE_TO_POLAR:
        remoteDir='/home/www/polar/global/class-1/archive/images/'
        subprocess.call('ssh tspindler@140.90.100.206 mkdir -p '+remoteDir+'/'+vDate,shell=True)
        subprocess.call('scp -q -r '+imageDir+'/* tspindler@140.90.100.206:'+remoteDir+'/'+vDate+'/.',shell=True)
    if CLEANUP_IMAGES:
        subprocess.call('rm -rf '+imageDir,shell=True)

