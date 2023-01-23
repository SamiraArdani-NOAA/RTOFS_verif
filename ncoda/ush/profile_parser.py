#!/bin/env python
"""
Profile Parser
version 1.0
Todd Spindler
30 August 2019

Notes
  30 Aug 2019 -- initial version in Git repo
  18 Dec 2020 -- added fixed buoys
  19 JUN 2022 -- ported to Hera (Gwen Chen)
  01 JUL 2022 -- changed the archive folders/links under /profiles on RZDM
                 to cut date instead of run date to be consistent with
                 the date in the Cut Date calendar on the webpage. Also
                 modified the orange text on top of the map. (Gwen Chen)  
"""

import os, sys
from datetime import datetime, timedelta
import pandas as pd
import subprocess

max_error = 4.0

outDir='/scratch2/NCEPDEV/stmp1/Lichuan.Chen/RTOFS-DA'
template="addMarker(ol.proj.fromLonLat([{},{}]),'id:\\t\\t{}\\nprofile:\\t{}\\nrcpt:\\t{}\\nerror:\\t{}','{}','profiles/{}/{}{:04n}{}');\n"

#--------------------------------------------------------------------
def write_js(data,fname,cutdate,imageType='.png'):
    if fname.lower().split('_')[0]=='buoy':
        imageType=fname.lower().split('_')[1]+'.gif'
    
    with open(outDir+'/markers.js','a') as f:
        f.write('function load_'+fname.lower()+'(){\nclearMarkers();\n')
        for val in data:
            f.write(template.format(val['lon'],
                                    val['lat'],
                                    val['id'],
                                    val['profile'],
                                    val['dtg'],
                                    val['stderr'],
                                    val['color'],
                                    cutdate.strftime('%Y%m%d'),
                                    fname.lower().split('_')[0],
                                    val['frame'],
                                    imageType))
        f.write("$('#group').html('"+fname.replace('_t',' Temperature').replace('_s',' Salinity')+"');\n")
        f.write("$('#cutdate').html(' ');\n")
        f.write("}\n")
        
#--------------------------------------------------------------------
def parse_section(lines,platform,pattern):
    data=[]
    section=[line.strip() for line in lines if (line.find('frame')>=0) & (line.find(pattern)>=0)]
    for line in section:
        prefix,id,line=line.split('"')
        frame=int(prefix.split()[1][:-1])
        profile=prefix.split()[3]
        
        if platform=='XBT':
            _,_,dtg,lat,lon,stderr=line.split()
        elif platform=='Argo':
            _,_,param,dtg,lat,lon,stderr=line.replace('BIAS','    ').split()
#            _,_,param,dtg,lat,lon,stderr=line.split()
        elif platform=='Glider':
            _,param,dtg,lat,lon,stderr=line.split()            
        elif platform=='TESAC':
            _,param,dtg,lat,lon,stderr=line.split()            
        elif platform=='Animal':
            _,_,param,dtg,lat,lon,stderr=line.split()            
        elif platform=='Buoy':
            _,_,param,dtg,lat,lon,stderr=line.split()
            
        stderr=float(stderr)
        if stderr < 3:
            color='blue'
        elif 3 <= stderr < 4:
            color='yellow'
        else: 
            color='red'

        if lon[-1]=='W':
            lon=-float(lon[:-1])
        else:
            lon=float(lon[:-1])
            
        if lat[-1]=='S':
            lat=-float(lat[:-1])
        else:
            lat=float(lat[:-1])
            
        data.append({'id':id.strip(),
                     'profile':profile,
                     'frame':frame,
                     'dtg':dtg,
                     'lon':str(lon),
                     'lat':str(lat),
                     'stderr':str(stderr),
                     'color':color
                    })
        
    return data

#--------------------------------------------------------------------
if __name__=='__main__':
    
    infile=sys.argv[1]    

    if os.path.exists('markers.js'):
        os.remove('markers.js')
    
    # read the out file
    with open(infile) as f:
        lines=f.readlines()
        
    cutdate = [x.strip() for x in lines if x.find('QC Cut Date Time Group:')>=0][0].split()[5]
    cutdate = datetime.strptime(cutdate,'%Y%m%d%H')
    
    # parse the sections
    xbt=parse_section(lines,'XBT','eXpendable BT')
    argo_t=parse_section(lines,'Argo','Argo Float Temp')
    argo_s=parse_section(lines,'Argo','Argo Float Salt')
    glider_t=parse_section(lines,'Glider','Glider Temp')
    glider_s=parse_section(lines,'Glider','Glider Salt')
    tesac_t=parse_section(lines,'TESAC','TESAC Temperature')
    tesac_s=parse_section(lines,'TESAC','TESAC Salinity')
    animal_t=parse_section(lines,'Animal','Animal Borne Temp')
    animal_s=parse_section(lines,'Animal','Animal Borne Salt')
    buoy_t=parse_section(lines,'Buoy','Fixed BUOY Temp')
    buoy_s=parse_section(lines,'Buoy','Fixed BUOY Salt')
    
    # extract buoy sublists to create gif loops
    buoy_t=pd.DataFrame.from_dict(buoy_t)
    buoy_s=pd.DataFrame.from_dict(buoy_s)
        
    for id in buoy_t.id.unique():
        loop=[f'{outDir}/png/buoy{int(x):04n}.png' for x in buoy_t[buoy_t.id==id].frame]
        loop = f"convert -delay 100 {' '.join(loop)} {outDir}/png/buoy{id}t.gif"
        subprocess.call(loop, shell=True)
        
    for id in buoy_s.id.unique():
        loop=[f'{outDir}/png/buoy{int(x):04n}.png' for x in buoy_s[buoy_s.id==id].frame]
        loop = f"convert -delay 100 {' '.join(loop)} {outDir}/png/buoy{id}s.gif"
        subprocess.call(loop, shell=True)
        
    buoy_t['stderr']=buoy_t.stderr.astype(float)
    buoy_s['stderr']=buoy_s.stderr.astype(float)
    maxerr_t=buoy_t.groupby('id').stderr.max()
    maxerr_s=buoy_s.groupby('id').stderr.max()
    buoy_t=buoy_t.drop_duplicates('id')
    buoy_s=buoy_s.drop_duplicates('id')
    buoy_t['stderr']=maxerr_t.values
    buoy_s['stderr']=maxerr_s.values
    buoy_t['frame']=buoy_t.id.astype(int)
    buoy_s['frame']=buoy_s.id.astype(int)
        
    buoy_t=[x._asdict() for x in buoy_t.itertuples(index=False)]
    buoy_s=[x._asdict() for x in buoy_s.itertuples(index=False)]
    
    datasets={'XBT':xbt,'Argo_t':argo_t,'Argo_s':argo_s,
              'Glider_t':glider_t,'Glider_s':glider_s,
              'TESAC_t':tesac_t,'TESAC_s':tesac_s,
              'Animal_t':animal_t,'Animal_s':animal_s,
              'Buoy_t':buoy_t,'Buoy_s':buoy_s}

    for key,data in datasets.items():
        write_js(data,key,cutdate)

