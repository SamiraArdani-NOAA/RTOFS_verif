"""
MMAB-specific classes and methods

author: Todd Spindler
version: 1.0

Notes
4/1/2019 - Corrected gridlines longitude misses edge case
3/1/2022 - Ported to Chen's Hera account
"""

import mpl_toolkits.basemap
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.image as image
import numpy as np
from datetime import datetime

class Basemap(mpl_toolkits.basemap.Basemap):
     #
     # Extend standard Basemap class with a new gridlines() method
     #
     def __init__(self,**kwargs):
          mpl_toolkits.basemap.Basemap.__init__(self,**kwargs)
          
     def gridlines(self,draw_labels=False,**kwargs):
          #
          #Draw nice gridlines on a map.
          #          
          nice=np.array((1,2,5,10,20,25,50,100))
          i=np.argmin(np.abs(nice-int(self.lonmax-self.lonmin)/8))
          j=np.argmin(np.abs(nice-int(self.latmax-self.latmin)/6))
          meridians=np.arange(0,361-nice[i],nice[i])
          parallels=np.arange(-90,91,nice[j])
          if self.projparams['proj']=='stere':
              self.drawparallels(parallels,labels=[1,1,0,0],**kwargs)
              self.drawmeridians(meridians,labels=[0,0,1,1],**kwargs)
          else:
              self.drawparallels(parallels,labels=[1,0,0,0],**kwargs)
              self.drawmeridians(meridians,labels=[0,0,0,1],**kwargs)
          if draw_labels:
               plt.xlabel('Longitude',**kwargs)
               plt.ylabel('Latitude',**kwargs)

class MidpointNormalize(colors.Normalize):
     def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
          self.midpoint = midpoint
          colors.Normalize.__init__(self, vmin, vmax, clip)
          
     def __call__(self, value, clip=None):
          result, is_scalar = self.process_value(value)
          (vmin,), _ = self.process_value(self.vmin)
          (vmax,), _ = self.process_value(self.vmax)
          resdat = np.asarray(result.data)
          result = np.ma.array(resdat, mask=result.mask, copy=False)
          x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
          res = np.interp(result, x, y)
          result = np.ma.array(res, mask=result.mask, copy=False)
          if is_scalar:
               result = result[0]
          return result

# set the colormap and centre the colorbar
class OldMidpointNormalize(colors.Normalize):
     """
     Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)
     
     e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
     """
     def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
          self.midpoint = midpoint
          colors.Normalize.__init__(self, vmin, vmax, clip)
        
     def __call__(self, value, clip=None):
          # I'm ignoring masked values and all kinds of edge cases to make a
          # simple example...
          x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
          return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))

def add_mmab_logos(alpha=1.0):
     """
     -----------------------------------------------
     add NOAA and NWS logos to an existing figure
     as well as some branding and dates
     -----------------------------------------------
     """
     noaa_logo=image.imread('/scratch1/NCEPDEV/stmp2/Samira.Ardani/Logos/NOAA_logo.png')
     nws_logo=image.imread('/scratch1/NCEPDEV/stmp2/Samira.Ardani/Logos/NWS_logo.png')
     fig=plt.gcf()
     fig.figimage(noaa_logo,alpha=alpha,
                  yo=fig.get_figheight()*fig.dpi-noaa_logo.shape[0])
     fig.figimage(nws_logo,alpha=alpha,
                  xo=fig.get_figwidth()*fig.get_dpi()-nws_logo.shape[1],
                  yo=fig.get_figheight()*fig.get_dpi()-nws_logo.shape[0])
     plt.annotate('NCEP/EMC/Verification Post Processing Product Generation Branch',
                  xy=(0.01,0.01),xycoords='figure fraction',
                  horizontalalignment='left',fontsize='x-small')
     plt.annotate(f'{datetime.now():%d %b %Y} on Hera',
                  xy=(0.99,0.01),xycoords='figure fraction',
                  horizontalalignment='right',fontsize='x-small')
    
def add_mmab_logos2(ax=None,alpha=1.0):
     # more branding, this time resizable
     noaa_logo = plt.imread('/scratch1/NCEPDEV/stmp2/Samira.Ardani/Logos/NOAA_logo.png')
     nws_logo = plt.imread('/scratch1/NCEPDEV/stmp2/Samira.Ardani/Logos/NWS_logo.png')
     if ax==None:
          ax=plt.gca()
     size=0.10
     llx,lly=ax.get_position().corners()[0]
     ulx,uly=ax.get_position().corners()[1]
     lrx,lry=ax.get_position().corners()[2]
     urx,ury=ax.get_position().corners()[3]
     
     # boilerplate at the bottom of main fig
     ax.annotate('NCEP/EMC/Verification Post Processing Product Generation Branch',
                 xy=(llx,lly/4.),xycoords='figure fraction',
                 horizontalalignment='left',verticalalignment='bottom',
                 fontsize='x-small')
     ax.annotate(f'{datetime.now():%d %b %Y} on $Hera$',
                 xy=(lrx,lry/4.),xycoords='figure fraction',
                 horizontalalignment='right',verticalalignment='bottom',
                 fontsize='x-small')
     # NOAA on the upper left
     plt.axes([ulx,uly+0.0,size,size])
     plt.axis('off')
     plt.imshow(noaa_logo,interpolation='gaussian',alpha=alpha)
     # NWS on the upper right
     plt.axes([urx-size,ury+0.0,size,size])
     plt.axis('off')
     plt.imshow(nws_logo,interpolation='gaussian',alpha=alpha)
     plt.sca(ax)
    
#---------------------------------------------------------------
# n-dimensional RMS calc                                        
# still old-style for rms(x,y)                                  
# new way for rms(t)                                            
# this normalizes the standard deviation by the sample size     
#---------------------------------------------------------------
def ndrms(obs,fcst,axis=None):
     """
     bias,rms,normsdev,cc,si=ndrms(obs,fcst,dim);
     
     bias and rms along a dimension
     if no dim then calculate rms, bias as timeseries
     """
     bias=[]
     rms=[]
     normsdev=[]
     cc=[]
     si=[]

     # this expects a three-dimensional array (time,lats,lons)
     # add the time dimension if it's not there
     
     if obs.ndim==2:
          obs=obs[np.newaxis,]
     if fcst.ndim==2:
          fcst=fcst[np.newaxis,]
          
     if axis != None:
          print('starting ndrms on axis',axis)
          obsmean=obs.mean(axis)
          fcstmean=fcst.mean(axis)
          B = fcstmean - obsmean
          B2 = B**2
          D = fcst-obs
          S=D.std(axis)
          S2=S**2
          bias=B.squeeze().astype('float32')
          rms=np.sqrt(S2 + B2).squeeze().astype('float32')
     else:
          print('starting taylor ndrms on ts with',obs.shape[0],'obs')
          for i in range(obs.shape[0]):
               print('processing obs',i)            
               tobs=obs[i,:,:]
               tfcst=fcst[i,:,:] 
               diff=tfcst-tobs
               B = tobs.mean() - tfcst.mean()
               obssdev=tobs.std()
               fcstsdev=tfcst.std()
               nsdev=fcstsdev/obssdev
               n = tobs.compressed().size
               E2 = (1/n)*(1/(obssdev**2))*np.ma.sum((tfcst-tobs)**2)
               RMS=E2**0.5        
               CC=np.ma.corrcoef(tobs.ravel(),tfcst.ravel())[0,1]
               SI=100.0*(((diff**2).mean())**0.5 - B**2)/tobs.mean()
               # collect values into timeseries
               bias.append(B)
               rms.append(RMS)
               cc.append(CC)
               normsdev.append(nsdev)
               si.append(SI)
     """
     else:
          print('starting ndrms on ts with',obs.shape[0],'obs')
          for i in range(obs.shape[0]):
               print('processing obs',i)            
               tobs=obs[i,:,:]
               tfcst=fcst[i,:,:]
               diff=tobs-tfcst
               obsmean=tobs.mean()
               fcstmean=tfcst.mean()
               B = fcstmean - obsmean
               B2 = B**2
               D = tfcst-tobs
               S=D.std()
               S2=S**2
               SI=100.0*(((diff**2).mean())**0.5 - B**2)/tobs.mean()
               bias.append(B)
               rms.append(np.sqrt(S2 + B2))
     """
     
     return bias,rms,normsdev,cc,si

#-------------------------------------------------------------------
# colorbar hack to resize it to the axis regardless of aspect ratios
#-------------------------------------------------------------------
def colorbar(mappable):
     ax = mappable.axes
     fig = ax.figure
     divider = make_axes_locatable(ax)
     cax = divider.append_axes("right", size="5%", pad=0.05)
     return fig.colorbar(mappable, cax=cax)
