"""
Florida Cable Transport Class-4 Validation System
Version 1.1  
Todd Spindler
5 June 2018

Changes
 6 Sep 2018 -- converted to Python3
               modified stats plots to draw blue nowcast line
24 Jul 2019 -- ported to Phase 3 (Mars)
 2 Aug 2020 -- return from Mars 
 1 Mar 2022 -- ported to Chen's Hera account with minor fixes
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as image
import matplotlib.dates as mdates
from matplotlib.dates import MonthLocator, DateFormatter
from matplotlib.dates import WeekdayLocator, MO, TU, WE, TH, FR, SA, SU
from netCDF4 import Dataset
import numpy as np
from scipy.integrate import trapz
from sklearn.metrics import mean_squared_error
from pyproj import Geod
import math
from datetime import datetime, timedelta
import pandas as pd
#import wget
#from ftplib import FTP
import multiprocessing as mp
import subprocess
import sys, os
import sqlite3
import warnings
#import ipdb; ipdb.set_trace()

# job control settings
WANT_POOL=True
ncores=12
UPDATE_DB=True
UPLOAD_TO_POLAR=False   ## if running on compute farms (cron), this wont work

# global subdirectories (will be env vars in next release)
baseDir='/scratch2/NCEPDEV/stmp1/Samira.Ardani/images/class-4/cable'
srcDir='/scratch2/NCEPDEV/ovp/Samira.Ardani'
dbfile=srcDir+'/VPPPG/Global_RTOFS/EMC_ocean-verification/cable/fix/cable_transports.db'
logoDir='/scratch2/NCEPDEV/ovp/Lichuan.Chen/Logos'
refDir=srcDir+'/VPPPG/Global_RTOFS/EMC_ocean-verification/cable/fix'
DCOMDir='/scratch2/NCEPDEV/ovp/Samira.Ardani/DCOM'
archDir='/scratch2/NCEPDEV/ovp/Samira.Ardani/Global/archive'
polarDir='/home/www/polar/global/class-4/cable/archive/images/'

#-----------------------------------------
# initialize the database if necessary    
#-----------------------------------------
def init_db(dbfile):

    table='cable'
    
    # connect to db file
    conn = sqlite3.connect(dbfile,detect_types=sqlite3.PARSE_DECLTYPES,timeout=30.0)
    c = conn.cursor()    
    
    # Create table if needed
    c.execute(
        'CREATE TABLE IF NOT EXISTS '+table\
            + ' (date timestamp,'\
            + ' forecast int,'\
            + ' bias real,'\
            + ' rmse real,'\
            + ' corr real,'\
            + ' scatter_index real,'\
            + ' unique(date,forecast))')
                                  
    # Save (commit) the changes
    conn.commit()
                 
    # close the connection
    conn.close()
    return
#-----------------------------------------------
# update the database                           
#-----------------------------------------------
def update_db(thedate,fcst,bias,rmse,corr,scatter_index):
    
    table='cable'
    
    init_db(dbfile)

    # connect to db file
    conn = sqlite3.connect(dbfile,detect_types=sqlite3.PARSE_DECLTYPES,timeout=30.0)
    c = conn.cursor()

    # Larger example that inserts many records at a time
    update = (thedate,fcst,bias,rmse,corr,scatter_index)
    
    c.execute(' REPLACE INTO '+table+' VALUES (?,?,?,?,?,?)', update)
    
    # Save (commit) the changes
    conn.commit()
                 
    # close the connection
    conn.close()    
    return
#-----------------------------------------------
# read from sqlite3 database into DataFrame     
#-----------------------------------------------
def read_db():
    init_db(dbfile)  # just in case it doesn't exist
    table='cable'
    conn = sqlite3.connect(dbfile,detect_types=sqlite3.PARSE_DECLTYPES,timeout=30.0)
    df=pd.read_sql_query(' SELECT DISTINCT date, forecast, bias, rmse, corr, scatter_index '+
                         ' FROM '+table+
                         ' ORDER BY date AND forecast',conn)
    return df
#-----------------------------------------------
# add logos to figure                           
#-----------------------------------------------
def mmab_logos(fig,alpha=1.0):
    # add some branding and dates
    noaa_logo=image.imread(logoDir+'/NOAA_logo.png')
    nws_logo=image.imread(logoDir+'/NWS_logo.png')
    fig.figimage(noaa_logo,alpha=alpha,
        yo=fig.get_figheight()*fig.dpi-noaa_logo.shape[0])
    fig.figimage(nws_logo,alpha=alpha,
        xo=fig.get_figwidth()*fig.get_dpi()-nws_logo.shape[1],
        yo=fig.get_figheight()*fig.get_dpi()-nws_logo.shape[0])
    plt.annotate('NCEP/EMC/Verification Post Processing Product Generation Branch',
    xy=(0.01,0.01),xycoords='figure fraction',
    horizontalalignment='left',fontsize='x-small')
    plt.annotate(f'{datetime.now():%d %b %Y} on $Hera$',xy=(0.99,0.01),xycoords='figure fraction',
    horizontalalignment='right',fontsize='x-small')

#-----------------------------------------------
# get cable transport data from AOML            
#-----------------------------------------------
def get_aoml(year=datetime.now().year):
    # get the AOML dataset for this year
    filename='FC_cable_transport_'+str(year)+'.dat'
    
    #url='ftp://ftp.aoml.noaa.gov/phod/pub/cable/'+filename
    #url='ftp://192.111.123.132/phod/pub/cable/'+filename

    #host='ftp.aoml.noaa.gov'
    host='192.111.123.132'
    remoteFile='/phod/pub/cable/'+filename
    
    # test the FTP server first then download the file
    """
    try:
        ftp=FTP(host,timeout=120)
        ftp.login()
        fn=wget.download('ftp://'+host+remoteFile,out=refDir+'/'+filename)
        # wget won't overwrite a file, so rename as needed
        if fn != refDir+'/'+filename:
            os.rename(fn,refDir+'/'+filename)                
    except:
        print('cannot FTP data from AOML -- continuing with current copy')
    """
    
    # read the AOML dataset
    names=['year','month','day','transport']
    cable=pd.read_csv(DCOMDir+'/AOML/'+filename,comment='%',names=names,delimiter=' ',
        skipinitialspace=True,header=None,usecols=list(range(4)))
    cable['date']=pd.to_datetime(cable[['year','month','day']])
    cable.index=cable.date
    cable['error']=2.0
    del cable['year'], cable['month'], cable['day'], cable['date']
    return cable
    
#-----------------------------------------------
# coordinate rotation                           
#-----------------------------------------------
def rotate(u,v,phi):
    # phi is in radians
    u2 =  u*math.cos(phi) + v*math.sin(phi)
    v2 = -u*math.sin(phi) + v*math.cos(phi)
    return u2,v2

#-----------------------------------------------
# full cross-section transport calculation      
#-----------------------------------------------
def calc_transport(dates,fcst):
    """
    Calculate the transport of water across the Florida Straits
    This extracts the section and integrates the flow through it.
    """
    transport=[]
    if fcst==0:
        fcst_str='n024'
    else:
        fcst_str='f{:03d}'.format(fcst)
    cable_loc=np.loadtxt(refDir+'/eightmilecable.dat',dtype='int',usecols=(0,1))
    eightmile_lat = 26.5167
    eightmile_lon = -78.7833%360
    wpb_lat = 26.7153425
    wpb_lon = -80.0533746%360
    cable_angle = math.atan((eightmile_lat-wpb_lat)/(eightmile_lon-wpb_lon))
    g=Geod(ellps='WGS84')
    
    for date in dates:
        print('processing',date.strftime('%Y%m%d'),'fcst',fcst)
        rundate=date-timedelta(fcst/24.)  # calc rundate from fcst and date
        ufile=archDir+'/'+rundate.strftime('%Y%m%d')+'/rtofs_glo_3dz_'+fcst_str+'_daily_3zuio.nc'
        vfile=archDir+'/'+rundate.strftime('%Y%m%d')+'/rtofs_glo_3dz_'+fcst_str+'_daily_3zvio.nc'
        
        try:
            udata=Dataset(ufile)
            vdata=Dataset(vfile)
        except:
            print(rundate,fcst,'not found -- continuing')
            transport.append(np.nan)
            continue
        
        lon=udata['Longitude'][:]
        lat=udata['Latitude'][:]
        depth=udata['Depth'][:]
                                        
        #ssh.set_fill_value(value=0.0)

        usection=np.zeros((depth.shape[0],cable_loc.shape[0]))
        vsection=np.zeros((depth.shape[0],cable_loc.shape[0]))
        
        #rows=np.array(cable_loc[:,0], dtype=np.intp)
        #cols=np.array(cable_loc[:,1], dtype=np.intp)
        #levs=np.arange(33, dtype=np.intp)
        
        #x[np.ix_(rows, columns)]
        
        udata=udata['u'][:].squeeze()
        vdata=vdata['v'][:].squeeze()
        
        for ncol,(row,col) in enumerate(cable_loc):
            usection[:,ncol]=udata[:,row,col].filled(fill_value=0.0)
            vsection[:,ncol]=vdata[:,row,col].filled(fill_value=0.0)
            
        lon=lon[cable_loc[:,0],cable_loc[:,1]]
        lat=lat[cable_loc[:,0],cable_loc[:,1]]

        # compute the distances along the track
        _,_,dist=g.inv(lon[0:-1],lat[0:-1],lon[1:],lat[1:])
        depth=np.diff(depth)
        usection=usection[:-1,:-1]
        vsection=vsection[:-1,:-1]

        dist,depth=np.meshgrid(dist,depth)        
        u,v=rotate(usection,vsection,cable_angle)        
        trans1=(v*dist*depth).sum()/1e6        
        #print(date.strftime('%Y-%m-%d'),' transport:',transport,'Sv')
        transport.append(trans1)
    
    return transport

#-----------------------------------------------
# retrieve model data                           
#-----------------------------------------------
def get_model(dates,fcsts):
    
    transport={'dates':dates}

    if WANT_POOL:
        pool=mp.Pool(processes=ncores)
        multiple_results=[]
        
    for fcst in fcsts:
        if WANT_POOL:
            multiple_results.append(pool.apply_async(calc_transport,(dates,fcst)))
        else:
            transport[fcst]=calc_transport(dates,fcst)
            
    if WANT_POOL:
        pool.close()
        pool.join()
        
        # process the results
        for fcst,res in zip(fcsts,multiple_results):
            #print(fcst)
            transport[fcst]=res.get(timeout=1)
    
    model=pd.DataFrame(transport)
    #model['validDates']=model.dates+timedelta(fcst/24.)
    #model.index=model.validDates
    model.index=model.dates
    del model['dates']
    #del model['validDates']

    return model

#-----------------------------------------------
# plot all transports                           
#-----------------------------------------------
def plot_transport(model,cable,fcsts,fcst,imageDir):
    """
    transport plotting routine
    """
    
    warnings.filterwarnings('ignore',module='matplotlib')

    # first off, calc stats and update the database as needed
    cable_bias,cable_rmse,cable_corr,cable_scatter_index=cable_stats(model,cable,fcst)
    
    # the date for this is the lesser of the cable data and model data
    if UPDATE_DB:
        update_db(min(model.index.max(),cable.index.max()).to_pydatetime(),fcst,cable_bias,cable_rmse,cable_corr,cable_scatter_index)
        
    if fcst==fcsts[0]:
        return
    
    # plotting 
    aoml_cable_error=2.0  # specified by AOML
    colors=['red','orange','yellow','green','blue','indigo','violet','brown']
    fig=plt.figure(dpi=150)
    ax=plt.subplot(111)
    cable.transport.plot(x_compat=True,color='black',xlim=datetime(2018,4,9),label='Cable',fontsize='x-small')
    plt.fill_between(cable.index,cable.transport+aoml_cable_error,cable.transport-aoml_cable_error,color='lightgrey',alpha=0.5,label='Cable Error')  
    #(cable.transport+aoml_cable_error).plot(x_compat=True,linestyle=':',color='black',label='Cable Error',fontsize='x-small')
    #(cable.transport-aoml_cable_error).plot(x_compat=True,linestyle=':',color='black',label='_none',fontsize='x-small')
    
    model[0].plot(x_compat=True,color='blue',linestyle='--',label='nowcast',fontsize='x-small')  # nowcast
    for ndx,fcst1 in enumerate(fcsts[1:]):
        if ndx==0:
            model[fcst1].plot(x_compat=True,linewidth=1,color='green',alpha=0.5,fontsize='x-small',label='all forecasts')
        else:
            model[fcst1].plot(x_compat=True,linewidth=1,color='green',alpha=0.5,fontsize='x-small',label='')
    model[fcst].plot(x_compat=True,linewidth=1,color='red',alpha=1.0,fontsize='x-small',label="{:03n} hr fcst".format(fcst))
    
    plt.xlim(xmin=model.index.min(),xmax=model.index.max())    
    plt.text(0,1,
             '{:03n} hr Fcst 7-day Stats\n'.format(fcst)+
             '         Bias = {:4.2f}\n'.format(cable_bias)+
             '    RMS Error = {:4.2f}\n'.format(cable_rmse)+
             '  Correlation = {:4.2f}\n'.format(cable_corr)+
             'Scatter Index = {:4.1f}'.format(cable_scatter_index),
             bbox=dict(facecolor='white', alpha=1.0),
             horizontalalignment='left',
             verticalalignment='top',
             fontsize='x-small',transform=ax.transAxes)
        
    #set ticks every week
    #ax.xaxis.set_major_locator(mdates.WeekdayLocator())
    #set major ticks format
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))    
    #fig.autofmt_xdate()
    plt.ylabel('Transport (Sv)',fontsize='x-small')
    plt.xlabel('Date',fontsize='x-small')
    plt.legend(fontsize='x-small',loc='lower left')
    plt.grid(which='both')
    plt.title('Florida Cable Transport (Sv)',fontsize='small')
    mmab_logos(fig)
    imageFile="cable_transport_{:03n}.png".format(fcst)
    plt.savefig(imageDir+'/'+imageFile)
    plt.close()
    return

#-----------------------------------------------
# compute 7-day stats                           
#-----------------------------------------------
def cable_stats(model,cable,fcst):
    
    # don't need the cable error
    
    if 'error' in list(cable.keys()):
        del cable['error']
    
    # combine only where both data sets are present
    both=pd.merge(cable,model,left_index=True,right_index=True,how='inner')
    both=both[both.index.max()-timedelta(7):]
    
    # not sure if this is a good idea or not
    both.dropna(inplace=True)
    
    diff=both[fcst] - both.transport
    bias=diff.mean()
    rmse=mean_squared_error(both.transport,both[fcst])**0.5
    if both[fcst].mean() != 0.0: 
        scatter_index=100.0*(((diff**2).mean())**0.5 - bias**2)/both.transport.mean() 
    else:
        scatter_index=np.nan
        
    corr=both[fcst].corr(both.transport)
    
    return bias, rmse, corr, scatter_index

#-----------------------------------------------
# plot 7-day stats                              
#-----------------------------------------------
def plot_stats(dates,imageDir):
    
    # set up date formats for plotting timeseries with once-per-day ticks
    loc=mdates.DayLocator()
    dateFmt=DateFormatter("%d %b %Y")

    # read in the database
    df=read_db()
    
    for param in ['bias','rmse','corr','scatter_index']:
        for ndx,fcst in enumerate(df.forecast.unique()):
            fig=plt.figure(dpi=150)
            ax=fig.add_subplot(111)
            df2=df[df.forecast==fcst].copy()
            #df2.index=df2.date
            for fcst2 in range(24,193,24):
                if fcst2==24:                    
                    df[df.forecast==fcst2].plot(ax=ax,x='date',y=param,x_compat=True,color='green',alpha=0.5,label='all forecasts',fontsize='x-small')
                else:
                    df[df.forecast==fcst2].plot(ax=ax,x='date',y=param,x_compat=True,color='green',alpha=0.5,label='',fontsize='x-small')
            df[df.forecast==0].plot(ax=ax,x='date',y=param,x_compat=True,color='blue',linestyle='--',label='nowcast',fontsize='x-small')
            
            df2.plot(ax=ax,x='date',y=param,x_compat=True,color='red',label="{:03n} hr fcst".format(fcst),fontsize='x-small')
            plt.xlim(xmin=dates.min(),xmax=dates.max())
            #plt.ylim(df[param].min(),df[param].max())
            plt.legend(fontsize='small')
            ax.xaxis.set_major_formatter(dateFmt)
            ax.xaxis.set_major_locator(loc)
            ax.set_xlabel('')
            plt.grid()
            plt.title("{} hr {} computed daily over previous week's data".format(fcst,param),fontsize='small')
            mmab_logos(fig)
            imageFile="cable_{}_{:03n}.png".format(param,fcst)
            plt.savefig(imageDir+'/'+imageFile)
            plt.close()
        
    return
#-------------------------------------------------------------------------
def rms(both,fcst):
    fcst=both[fcst]
    obs=both['transport']
    B = fcst.mean() - obs.mean()
    B2 = B**2
    D = fcst-obs
    S=D.std()
    S2=S**2
    return np.sqrt(S2 + B2)

#-------------------------------------------------------------------------
def bias(both,fcst):
    fcst=both[fcst]
    obs=both['transport']
    return fcst.mean() - obs.mean()

#-------------------------------------------------------------------------
if __name__ == "__main__":
    
    DateSet=False
    if len(sys.argv)>1:
        want_date=datetime.strptime(sys.argv[1],'%Y%m%d')
        DateSet=True
        
    print('Starting Cable V&V at',datetime.now())
        
    fcsts=list(range(0,193,24))
        
    # get the latest cable data first and reference the model data from that
    cable_lastyear=get_aoml(datetime.now().year-1)
    cable_now=get_aoml()
    cable=cable_lastyear.append(cable_now)
    #cable=cable_now.copy()
    
    # read in the database

    df=read_db()
    
    # decide if to bring db up to date or backfill a single date
    if not DateSet:
        stop_date=cable.index.max()
        if df.size > 0:            
            start_date=df.date.max()
        else:
            start_date=stop_date
    else:
        start_date=want_date
        stop_date=want_date
        cable=cable[:stop_date] # cut off cable at stop date
        
    # loop through from last db date to last cable date
    for end_date in pd.date_range(start_date,stop_date):
        
        dates=pd.date_range(end=end_date,periods=60)
        model=get_model(dates,fcsts)
        imageDir=baseDir+'/'+end_date.strftime('%Y%m%d')
        if not os.path.isdir(imageDir):
            os.makedirs(imageDir)        

        for fcst in fcsts:
            plot_transport(model,cable,fcsts,fcst,imageDir)
            
        plot_stats(dates,imageDir)

        if UPLOAD_TO_POLAR:
            remoteDir=polarDir+end_date.strftime('%Y%m%d')
            subprocess.call('ssh emc.rtofs@140.90.100.206 mkdir -p '+remoteDir,shell=True)
            subprocess.call('scp '+imageDir+'/*.png emc.rtofs@140.90.100.206:'+remoteDir,shell=True)
            subprocess.call('rm -rf '+imageDir,shell=True)
    
