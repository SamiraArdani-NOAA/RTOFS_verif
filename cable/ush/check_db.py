"""
 walkthrough cable_transport.py manually
"""
import cable_transport as cb
from datetime import datetime
import pandas as pd
import os

baseDir='/scratch2/NCEPDEV/stmp1/Deanna.Spindler/images/class-4/cable'

want_date=datetime.strptime('20220105','%Y%m%d')
fcsts=list(range(0,25,24))

cable_lastyear=cb.get_aoml(datetime.now().year-1)
cable_now=cb.get_aoml()
cable=cable_lastyear.append(cable_now)

df=cb.read_db()

start_date=want_date
stop_date=want_date
cable=cable[:stop_date]

for end_date in pd.date_range(start_date,stop_date):
	dates=pd.date_range(end=end_date,periods=60)
	model=cb.get_model(dates,fcsts)
	imageDir=baseDir+'/'+end_date.strftime('%Y%m%d')
	if not os.path.isdir(imageDir):
		os.makedirs(imageDir)

	#cb.plot_transport(model,cable,fcsts,24,imageDir)

	
	cable_bias,cable_rmse,cable_corr,cable_scatter_index=cb.cable_stats(model,cable,fcst)

	
all_transports=cable.merge(model,left_index=True,right_index=True)
# "transport" column is from the AOML cable
# columns with 0, 24, 48... are transports from the model fcst
dbfile='all_cable_transports.db'
conn = sqlite3.connect(dbfile,detect_types=sqlite3.PARSE_DECLTYPES,timeout=30.0)
all_transports.to_sql('all_transports',conn,if_exists='append')
conn.close()

df2=df.copy()
df2.rename(columns={'date':'dates'},inplace=True)
df2.set_index('dates',inplace=True)

df2.merge(model


df2['model_T']=NaN
df2[df2.forecast==0].model_T=model[0]

df2.loc[df2.index==model.index,df2.forecast==0].model_T=model[0]

