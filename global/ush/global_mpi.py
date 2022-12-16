#!/bin/env python
"""
title: Global RTOFS Graphics Visualization with MPI
author: Todd Spindler
version: 1.3

Notes: 

This version maps T,S anomaly fields in addition to all other params

10 Jul 2018 (v1.1): Converted to Python 3. Changed regions diction to 
                    OrderedDict to fix wandering color limits.  PIL image
                    palette changed to WEB, reducing file sizes.
14 Jun 2019 (v1.2): removed OrderedDict as all dictionaries in Python3.7 are
		    ordered by default.  Moved all code to vlab/git and
                    VPPPG_Marine_Dev directory.  Modified climo load step
                    (load_woa) to interpolate to day.
22 Jul 2019 (v1.3): ported to Phase 3 (Mars)
 7 Aug 2020       : ported to Hera
15 Dec 2020 (v1.4): ported to RTOFSv2
 1 Mar 2022       : ported to Chen's Hera account with minor fixes
"""

from ipyparallel import Client
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import pandas as pd
from datetime import date, datetime, timedelta
import time
import sys, os
import pickle
import subprocess
import io

wantMPI=True

# set numpy error handling for overflow conditions
np.seterr(over='ignore')
    
#------------------------------------------------------------------------
def loadData(param,ncDate,fcst):

#    ncDir='/scratch2/NCEPDEV/ovp/Lichuan.Chen/Global/archive/'
    ncDir='/scratch2/NCEPDEV/ovp/Lichuan.Chen/Global/archive/'

    nc={           'temperature':'rtofs_glo_3dz_'+fcst+'_daily_3ztio.nc',
                      'salinity':'rtofs_glo_3dz_'+fcst+'_daily_3zsio.nc',
                             'u':'rtofs_glo_3dz_'+fcst+'_daily_3zuio.nc',
                             'v':'rtofs_glo_3dz_'+fcst+'_daily_3zvio.nc',
                             'ssh':'rtofs_glo_2ds_'+fcst+'_diag.nc',
         'mixed_layer_thickness':'rtofs_glo_2ds_'+fcst+'_diag.nc',
                  'ice_coverage':'rtofs_glo_2ds_'+fcst+'_ice.nc',
                 'ice_thickness':'rtofs_glo_2ds_'+fcst+'_ice.nc'}
    try:
        nco=Dataset(ncDir+ncDate+'/'+nc[param])
    except RuntimeError as e:
        print('*******************************')
        print('I/O error: {0}'.format(e)       )
        print('Stopping now'                   )
        print('*******************************')
        sys.exit()
                                                                
    lon=nco.variables['Longitude'][:]
    lat=nco.variables['Latitude'][:]
    data=nco.variables[param]
    lon[-1,]=lon[-2,]  # the usual patch
    return lon,lat,data

#------------------------------------------------------------------------
def drawMap(param,nick,region,lon,lat,data,u,v,lonlim,latlim,title,imageDir,imagefile,clim):
    import matplotlib as mpl
    mpl.use('Agg')
    import numpy as np
    import numpy.ma as ma
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from mmab_toolkit import Basemap, MidpointNormalize
    from PIL import Image
    import seaborn as sns
    import cmaps
    import io
    import sys
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
    #-------------------------------------
    def add_mmab_logos2(ax=None,alpha=1.0):
        # more branding, this time resizable
        noaa_logo = plt.imread('/scratch2/NCEPDEV/ovp/Lichuan.Chen/Logos/NOAA_logo.png')
        nws_logo = plt.imread('/scratch2/NCEPDEV/ovp/Lichuan.Chen/Logos/NWS_logo.png')
        size=0.10
        if ax==None:
            ax=plt.gca()        
        # boilerplate at the bottom of main fig
        ax.annotate('NCEP/EMC/Verification Post Processing Product Generation Branch',
                    xy=(0,0),xycoords='figure fraction',
                    horizontalalignment='left',verticalalignment='bottom',
                    fontsize='x-small')
        ax.annotate(f'{datetime.now():%d %b %Y} $on Hera$',
                    xy=(1,0),xycoords='figure fraction',
                    horizontalalignment='right',verticalalignment='bottom',
                    fontsize='x-small')
        # NOAA on the upper left
        #ax2=plt.axes([ulx,uly+size/2.,size,size],anchor='NW')
        ax2=plt.axes([0,1-size,size,size],anchor='NW')
        plt.imshow(noaa_logo,interpolation='gaussian',alpha=alpha)
        plt.axis('off')
        # NWS on the upper right
        #ax3=plt.axes([urx-size,ury+size/2.,size,size],anchor='NE')
        ax3=plt.axes([1-size,1-size,size,size],anchor='NE')
        plt.imshow(nws_logo,interpolation='gaussian',alpha=alpha)
        plt.axis('off')        
        plt.sca(ax)
        return
    
    #------------------------------
    plt.style.use('classic')
    fig=plt.figure()
            
    if nick == 'indian':
        lon=np.concatenate((lon[:,-1000:]-360,lon[:,1:-1000]),1)
        lat=np.concatenate((lat[:,-1000:],lat[:,1:-1000]),1)
        data=ma.concatenate((data[:,-1000:],data[:,1:-1000]),1)
            
    if nick == 'southern':
        m = Basemap(projection='spstere', 
            lon_0=270, 
            lat_0=-90, 
            boundinglat=region['lat'][1], 
            round=True, 
            resolution=region['resol'])
            
    elif nick == 'arctic':
        m = Basemap(projection='npstere',
            lon_0=270,
            lat_0=90,
            boundinglat=region['lat'][0],
            round=True,
            resolution=region['resol'])

    else:            
        m=Basemap(projection='mill',
            llcrnrlon=lonlim[0],
            urcrnrlon=lonlim[1], 
            llcrnrlat=latlim[0],
            urcrnrlat=latlim[1],
            resolution=region['resol'])

    x,y=m(lon,lat)

    # set the plot limits as needed
    if param == 'salinity':
        vmin,vmax=30,clim[1]
    elif param == 'mixed_layer_thickness':
        vmin,vmax=clim[0],min(clim[1],500)
    elif param == 'curr' and nick=='global':
        vmin,vmax=clim[0],(clim[1]-clim[0])*2/3+clim[0]
    else:
        vmin,vmax=clim[0],clim[1]

    # create nice bicolor map for some params
    if param=='temp_anom' or param=='salin_anom' or param=='ssh':
        #cmap=ListedColormap(sns.color_palette("RdBu_r",256))
        cmap=cmaps.NCV_jaisnd
        norm=MidpointNormalize(midpoint=0.0,vmin=vmin,vmax=vmax)
    else:
        cmap=cmaps.WhViBlGrYeOrRe
        norm=None
        
    # now get out the crayolas
    cs=m.pcolormesh(x,y,data,shading='flat',cmap=cmap,vmin=vmin,vmax=vmax,norm=norm)
        
    #plt.clim(vmin,vmax)
        
    if param =='curr' and nick != 'global':
        rskip=int(np.floor(data.shape[0]/50))
        cskip=int(np.floor(data.shape[1]/50))
        m.quiver(x[0::rskip,0::cskip],y[0::rskip,0::cskip], \
                 u[0::rskip,0::cskip],v[0::rskip,0::cskip], \
                 color='white',pivot='middle',scale_units='xy', \
                 alpha=0.7)
             
    m.drawcoastlines()
    m.fillcontinents()
    m.gridlines()
    m.drawmapboundary(linewidth=2.0)
    add_mmab_logos2()
    
    if nick == 'eqpac' or nick == 'med' or nick == 'spac':
        csize='7%'
        cpad='10%'
    else:
        csize='5%'
        cpad='7%'
        
    if param == 'salinity':
        plt.clim(vmin,data.max())
        cbar=m.colorbar(cs,location='bottom',extend='min',size=csize,pad=cpad)
    elif param == 'mixed_layer_thickness':
        cbar=m.colorbar(cs,location='bottom',extend='max',size=csize,pad=cpad)
    elif param == 'curr' and nick=='global':
        cbar=m.colorbar(cs,location='bottom',extend='max',size=csize,pad=cpad)
    else:
        cbar=m.colorbar(cs,location='bottom',size=csize,pad=cpad)
    
    cbar.ax.tick_params(labelsize='x-small')
    plt.title(title,fontsize='small')

    saveImage(imageDir+'/large/'+imagefile,dpi=150)
    saveImage(imageDir+'/small/'+imagefile,dpi=80)    
    plt.close()
    
    return
#------------------------------------------------------------------------
def initialize():
    regions = {}    
    regions['global']={'name':'Global',
        'resol':'c','lat':(-80,90),'lon':(74,434),
        'xlim':(1,4500),'ylim':(1,3298)}    
    regions['natl'] = {'name':'North Atlantic',
        'resol':'c','lat':(0,60),'lon':(260,390),
        'xlim':(2324,4500),'ylim':(1505,3115)}
    regions['satl'] = {'name':'South Atlantic',
        'resol':'c','lat':(-80,0),'lon':(280,390),
        'xlim':(2574,3949),'ylim':(1,1505)}
    regions['eqatl'] = {'name':'Equatorial Atlantic',
        'resol':'c','lat':(-30,30),'lon':(280,390),
        'xlim':(2574,3949),'ylim':(1112,1898)}
    regions['npac'] = {'name':'North Pacific',
        'resol':'c','lat':(0,70),'lon':(95,290),
        'xlim':(262,4500),'ylim':(1505,3297)}
    regions['spac'] = {'name':'South Pacific',
        'resol':'c','lat':(-80,0),'lon':(95,295),
        'xlim':(262,2761),'ylim':(1,1505)}
    regions['eqpac'] = {'name':'Equatorial Pacific',
        'resol':'c','lat':(-30,30),'lon':(95,295),
        'xlim':(262,2761),'ylim':(1112,1898)}
    regions['southern'] = {'name':'Southern Ocean',
        'resol':'c','lat':(-80,-30),'lon':(74,434),
        'xlim':(1,4500),'ylim':(1,1249)}
    regions['indian'] = {'name':'Indian Ocean',
        'resol':'c','lat':(-75,30),'lon':(10,130),
        'xlim':(1,4500),'ylim':(1,1898)}    
    regions['arctic'] = {'name':'Arctic Ocean',
        'resol':'c','lat':(47,90),'lon':(74,434),
        'xlim':(1,4500),'ylim':(2173,3298)}    
    regions['med'] = {'name':'Mediterranean Sea',
        'resol':'l','lat':(29,48),'lon':(350,405),
        'xlim':(3449,4136),'ylim':(1885,2209)}
    regions['gulfstream'] = {'name':'Gulf Stream Region',
        'resol':'l','lat':(15,47),'lon':(260,306),
        'xlim':(2324,2899),'ylim':(1695,2172)}
    regions['kuroshio'] = {'name':'Kuroshio Region',
        'resol':'l','lat':(15,45),'lon':(115,155),
        'xlim':(512,1011),'ylim':(1695,2136)}
    regions['agulhas'] = {'name':'Agulhas Region',
        'resol':'l','lat':(-50,-15),'lon':(360,405),
        'xlim':(3574,4136),'ylim':(782,1315)}
    regions['somalia'] = {'name':'Somalia Region',
        'resol':'l','lat':(-35,15),'lon':(385,420),
        'xlim':(3887,4324),'ylim':(1038,1694)}
    regions['northbrazil'] = {'name':'North Brazil Region',
        'resol':'l','lat':(-15,15),'lon':(295,340),
        'xlim':(2762,3324),'ylim':(1316,1694)}
    regions['persiangulf'] = {'name':'Somalia & Persian Gulf',
        'resol':'l','lat':(0,35),'lon':(400,425),
        'xlim':(4074,4386),'ylim':(1505,1972)}
    regions['australia'] = {'name':'Australia & New Zealand',
        'resol':'l','lat':(-50,-10),'lon':(110,180),
        'xlim':(449,1324),'ylim':(782,1379)}
    regions['indonesia'] = {'name':'Indonesia & Philippines',
        'resol':'l','lat':(-10,25),'lon':(100,130),
        'xlim':(324,699),'ylim':(1380,1827)}
    regions['gulfmex'] = {'name':'Gulf of Mexico',
        'resol':'l','lat':(15,35),'lon':(260,285),
        'xlim':(2324,2636),'ylim':(1695,1972)}
    regions['westconus'] = {'name':'Western CONUS',
        'resol':'l','lat':(10,71),'lon':(202,282),
        'xlim':(1454,2881),'ylim':(1631,3298)}
    regions['alaska'] = {'name':'Alaska',
        'resol':'l','lat':(46,75),'lon':(180,218),
        'xlim':(1100,2500),'ylim':(2000,3500)}

    params     = ('temperature','salinity','ssh','mixed_layer_thickness','ice_coverage','ice_thickness','curr',   'temp_anom',   'salin_anom')
    paramNames = ('Temperature','Salinity','SSH','Mixed Layer Thickness','Ice Coverage','Ice Thickness','Current','Temp Anomaly','Salin Anomaly')
    paramUnits = ('deg C',      'psu',     'm',  'm',                    'pct',         'm',            'm/s',    'deg C',       'psu')
    
    #params     = ('temp_anom',)
    #paramNames = ('Temp Anomaly',)
    #paramUnits = ('deg C',)

    zindex=(0,  6,  9, 13, 15,  18,  25,  27,  29)
    #zindex=(0,)
    
    depths=(0,100,200,500,700,1000,2000,3000,4000)
    #depths=(0,)
        
    return regions,params,paramNames,paramUnits,zindex,depths
#------------------------------------------------------------------------
def startCluster(numEngines,fcst):
    
    print('starting cluster')
    subprocess.Popen(["ipcluster",
        "start",
        "--n=" + str(numEngines),
        "--profile=mpi",
        "--ip='*'",                      
        "--cluster-id="+fcst,
        "--daemonize"])
        
    time.sleep(60) # wait for the cluster to spin up

    print('initializing engines')
    rc=Client(profile='mpi',cluster_id=fcst)

    return rc
#------------------------------------------------------------------------
# subprocess version of cluster shutdown
#------------------------------------------------------------------------
def stopCluster(fcst):

    # shut down cluster
    print('Shutting down cluster')
    subprocess.Popen(["ipcluster", 
        "stop",
        "--profile=mpi",
        "--cluster-id="+fcst])
        
    time.sleep(60) # wait for the cluster to shut down
    print('Cluster should be stopped')

    return
#-------- subregion plot limits ----------------------------
def getDataLimits(clim,data):
    return [min(data.min(),clim[0]),max(data.max(),clim[1])]    
#------------------------------------------------
def get_woa(param,ncDate,fcst):
    """
    Read World Ocean Atlas 2013 climatology data        
    WOA monthly data is centered on mid-month, need to scale the model date 
    against the date range to create a weighted average of the straddling 
    month fields.
    """    
    climoDir='/scratch2/NCEPDEV/ovp/Lichuan.Chen/Global/climo/WOA/regrid'
    if fcst=='n024':
        fcst_sign=-1
    elif fcst=='n048':
        fcst='n000'
        fcst_sign=-1
    else:
        fcst_sign=1
    climoDate=datetime.strptime(ncDate,'%Y%m%d')+timedelta(fcst_sign*int(fcst[1:])/24.)        
    if climoDate.day==15:  # even for Feb, just because
        climofile="woa13_regrid_{:02n}.nc".format(climoDate.month)
        data=Dataset(climoDir+'/'+climofile)
        data=data[param][:]  # surface only
    else:  # need to scale things
        if climoDate.day < 15:
            start=climoDate-pd.Timedelta(days=20)
            start=pd.Timestamp(start.year,start.month,15)
            stop=pd.Timestamp(climoDate.year,climoDate.month,15)
        else:
            start=pd.Timestamp(climoDate.year,climoDate.month,15)
            stop=climoDate+pd.Timedelta(days=20)
            stop=pd.Timestamp(stop.year,stop.month,15)
        left=(climoDate-start)/(stop-start)
        #right=(stop-climoDate)/(stop-start)
        climofile1="woa13_regrid_{:02n}.nc".format(start.month)
        climofile2="woa13_regrid_{:02n}.nc".format(stop.month)
        data1=Dataset(climoDir+'/'+climofile1)
        data2=Dataset(climoDir+'/'+climofile2)
            
        data1=data1[param][:]
        data2=data2[param][:]
        data=data1+((data2-data1)*left)

    return data
#------------------------------------------------------------------------
if __name__ == '__main__':
    
    t = time.time()

    ncDate=sys.argv[1]
    fcst=sys.argv[2]
    numEngines=int(sys.argv[3])
    print("Arguments:",ncDate,fcst,numEngines)
    
    print('Starting Global MPI Visualization at',datetime.now(),'for',ncDate)
    
    if wantMPI:        
        rc=startCluster(numEngines,ncDate+'_'+fcst)
        # check that cluster started properly
        print('Cluster health check:')
        print('numEngines =',numEngines,'Actual Engines =',len(rc.ids))
        if len(rc.ids) != numEngines:
            print('WARNING: Cluster IDs do not match requested number')
            print('Restarting cluster')
            stopCluster(fcst)
            rc=startCluster(numEngines,ncDate+'_'+fcst)
            print('numEngines =',numEngines,'Actual Engines =',len(rc.ids))

        rc[:].use_cloudpickle()
        view=rc.load_balanced_view()        
        view.set_flags(retries=1)
        view.block=True
        
    # initialize regions, parameters and depths
    regions,params,paramNames,paramUnits,zindex,depths = initialize()

    # colorbar limits array
    new_clim=np.zeros([len(depths),len(regions),len(params),2])
    new_clim[:,:,:,0]=999.
    
    #unpickle the clim array if available
    try:
        f=open('/scratch2/NCEPDEV/ovp/Lichuan.Chen/VPPPG/Global_RTOFS/EMC_ocean-prod-gen/z-level/fix/clim.pkl','rb')
        old_clim=pickle.load(f)
        f.close()
    except:
        old_clim=new_clim
        
    print('starting on',ncDate,'forecast',fcst)
    result=[]
        
    # loop through parameters
    for nparam, (param,pname,punit) in enumerate(zip(params,paramNames,paramUnits)):
        print('starting on parameter',param)
        if param=='u' or param == 'v':
            continue
        #print('loading',param)
        if param == 'curr':
            lon,lat,udata=loadData('u',ncDate,fcst)
            lon,lat,vdata=loadData('v',ncDate,fcst)
        elif param == 'temp_anom' or param == 'salin_anom':
            if param == 'temp_anom':
                lon,lat,data=loadData('temperature',ncDate,fcst)
                climo=get_woa('temperature',ncDate,fcst)
            else:
                lon,lat,data=loadData('salinity',ncDate,fcst)
                climo=get_woa('salinity',ncDate,fcst)
        else:
            lon,lat,data=loadData(param,ncDate,fcst)
        
        # loop through levels
        for nlevel,(z,depth) in enumerate(zip(zindex,depths)):
            if depth > 0 and param != 'temperature' and \
                             param != 'salinity' and  \
                             param != 'curr' and  \
                             param != 'temp_anom' and \
                             param != 'salin_anom':
                break
            if depth > 1000 and \
                (param == 'temp_anom' or param == 'salin_anom'):
                break
            if param == 'curr':
                ulev=udata[:,z,].squeeze()
                vlev=vdata[:,z,].squeeze()
                levData=(ulev**2+vlev**2)**0.5
            else:
                if  param == 'temperature' or \
                    param == 'salinity' or \
                    param == 'curr' or \
                    param == 'temp_anom' or \
                    param == 'salin_anom':
                    levData=data[:,z,].squeeze()
                    if param=='temp_anom' or param=='salin_anom':
                        levData=levData - climo[nlevel,]
                else:
                    levData=data[:].squeeze()
                    
            # loop through regions
            for nreg,(nick,region) in enumerate(regions.items()):
            
                # skip ice in warm regions, mask zero-values otherwise
                if (param == 'ice_coverage' or param == 'ice_thickness') and \
                    (nick == 'med' or \
                    nick == 'eqatl' or \
                    nick == 'eqpac' or \
                    nick == 'gulfstream' or \
                    nick == 'northbrazil' or \
                    nick == 'gulfmex' or \
                    nick == 'kuroshio' or \
                    nick == 'australia' or \
                    nick == 'indonesia' or \
                    nick == 'persiangulf' or \
                    nick == 'agulhas' or \
                    nick == 'somalia'):
                    continue
                else:
                    levData.mask = ma.mask_or(levData.mask,levData==0.0)

                xlim=region['xlim']
                ylim=region['ylim']
                lons=region['lon']
                lats=region['lat']
                lon2=lon[ylim[0]:ylim[1],xlim[0]:xlim[1]]
                lat2=lat[ylim[0]:ylim[1],xlim[0]:xlim[1]]
                data2=levData[ylim[0]:ylim[1],xlim[0]:xlim[1]]
                
                if param == 'curr':
                    u2=ulev[ylim[0]:ylim[1],xlim[0]:xlim[1]]
                    v2=vlev[ylim[0]:ylim[1],xlim[0]:xlim[1]]
                    
                if param == 'ice_coverage' or param == 'ice_thickness':
                    data2.mask = ma.mask_or(data2.mask,data2==0.0)

                # data for f192 currents is very strange at the moment
                dataMin=str(data2.min())
                dataMax=str(data2.max())
                if dataMin!='--':
                    dataMin='{:6.2f}'.format(data2.min())
                if dataMax!='--':
                    dataMax='{:6.2f}'.format(data2.max())
                
                title=region['name']+' '+pname+' ('+punit+')\n'+ \
                    ncDate+' '+fcst+'\nDepth:'+'{: >5}'.format(depth)+' m'+ \
                    ' min: ' + dataMin + ' max: ' + dataMax
                    
                imageDir='/scratch2/NCEPDEV/stmp1/Lichuan.Chen/images/mpi/'+ncDate

                image='rtofs_'+nick+'_'+param+'_'+fcst+'_'+ \
                    '{0!s:0>3}'.format(depth)+'.png'

                # update plot limits from 000Z nowcast
                new_clim[nlevel,nreg,nparam,:]=getDataLimits(\
                    old_clim[nlevel,nreg,nparam,],data2)
                    
                if wantMPI:
                    if param == 'curr':
                        result.append(view.apply_async(drawMap,
                            param,nick,region,lon2,lat2,data2,
                            u2,v2,lons,lats,title,imageDir,image,
                            old_clim[nlevel,nreg,nparam,]))
                    else:
                        result.append(view.apply_async(drawMap,
                            param,nick,region,lon2,lat2,data2,
                            None,None,lons,lats,title,imageDir,image,
                            old_clim[nlevel,nreg,nparam,]))
                else:
                    if param =='curr':
                        drawMap(param,nick,region,lon2,lat2,data2,
                            u2,v2,lons,lats,title,imageDir,image,
                            old_clim[nlevel,nreg,nparam,])
                    else:
                        drawMap(param,nick,region,lon2,lat2,data2,
                            None,None,lons,lats,title,imageDir,image,
                            old_clim[nlevel,nreg,nparam,])
                        
    if wantMPI:
        print('waiting for results')
        view.wait(result)
        result=[]
        # shut down cluster
        subprocess.Popen(["ipcluster", 
                    "stop",
                    "--profile=mpi",
                    "--cluster-id="+ncDate+'_'+fcst])

    if fcst == 'n048':
        f=open('/scratch2/NCEPDEV/ovp/Lichuan.Chen/VPPPG/Global_RTOFS/EMC_ocean-prod-gen/z-level/fix/clim.pkl','wb')
        pickle.dump(new_clim,f)
        f.close()
                                                
    print('Elapsed: %.2f seconds'%(time.time() - t))
    
