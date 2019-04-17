import iris 
import iris.plot as iplt
import iris.quickplot as qplt
import iris.coord_categorisation
from iris.util import *

import iris.analysis.cartography 
import cartopy.crs as ccrs
import iris.coord_categorisation as iriscc
import cartopy.feature as cfe
from mpl_toolkits.basemap import Basemap, maskoceans

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import warnings
from shiftedColorMap import*
from iris.experimental.equalise_cubes import *

import glob
import os 
from itertools import groupby

import pandas as pd


#################################### COMPUTE NAO INDEX  ############################################################################################

def NAO(data_in, seasonal):


 Lisbon=[('latitude', 38.72 ), ('longitude',  -9.13)]
 Ryk= [('latitude', 64.14 ), ('longitude',  -21.94)]
 PDelg= [('latitude', 37.73 ), ('longitude',  -25.67)]
 

 if seasonal == 1:

  North=data_in.interpolate(Ryk, iris.analysis.Linear())
  South=data_in.interpolate(Lisbon, iris.analysis.Linear())

  #compute DJFM winter means
  months=data_in.shape[0]
  North_DJFM=iris.cube.CubeList() 
  South_DJFM=iris.cube.CubeList() 

  for i in range (12, months-12, 12):
     North_DJFM.append(North[i-1:i+3].collapsed('time', iris.analysis.MEAN))
     South_DJFM.append(South[i-1:i+3].collapsed('time', iris.analysis.MEAN))
  North= North_DJFM.merge_cube()
  South= South_DJFM.merge_cube()
  #print North_DJFM[0]
 
 else: #for annual means
  North=data_in.interpolate(Ryk, iris.analysis.Linear())
  South=data_in.interpolate(PDelg, iris.analysis.Linear())


 #long term mean to be used to calculate anomalies
 North_mean=North.collapsed('time', iris.analysis.MEAN)
 North_stdev=North.collapsed('time', iris.analysis.STD_DEV)

 South_mean=South.collapsed('time', iris.analysis.MEAN)
 South_stdev=South.collapsed('time', iris.analysis.STD_DEV)

 #normalize each month separately
 North_stnd=(North - North_mean)/North_stdev
 South_stnd=(South - South_mean)/South_stdev

 NAO= South_stnd - North_stnd

 return NAO


################################################## PLOt NAO INDEX  ###########################################################################################

def Plot(NAO_MO, NAO_AR, title_1, title_2):

 mean_AR=NAO_AR.collapsed('time', iris.analysis.MEAN)
 stdev_AR=NAO_AR.collapsed('time', iris.analysis.STD_DEV)
 mean_MO=NAO_MO.collapsed('time', iris.analysis.MEAN)
 stdev_MO=NAO_MO.collapsed('time', iris.analysis.STD_DEV)
 print mean_AR.data, mean_MO.data, stdev_AR.data, stdev_MO.data

 months=NAO_AR.shape[0]
 #times= pd.date_range(start='1850-01-01', periods=months, freq='AS') 
 times=np.arange(0,months,1)

 zero_line=np.empty(months); zero_line.fill(0)

 fig = plt.figure()

 ax1 = fig.add_subplot(211)

 fig.autofmt_xdate()
 plt.plot(times, zero_line, c='black' )
 plt.plot(times, NAO_AR.data, c='black',linewidth=1, label='ARCHER')
 plt.fill_between(times, 0, np.ma.masked_where(NAO_AR.data <= 0, NAO_AR.data) , alpha=0.5, facecolor='red')
 plt.fill_between(times, 0, np.ma.masked_where(NAO_AR.data >= 0, NAO_AR.data) , alpha=0.5, facecolor='blue')
 plt.xlabel('Time(years)', fontsize=14)
 plt.ylabel('NAO Index', fontsize=14)
 plt.title(title_1, fontsize=16)
 #plt.ticklabel_format(axis='y',useOffset=False)  #only on y axis if we are using a date array created with pandas
 #ax.xaxis.set_major_locator(mdates.YearLocator(20))
 plt.legend()


 ax2 = fig.add_subplot(212)

 fig.autofmt_xdate()
 plt.plot(times, zero_line, c='black' )
 plt.plot(times, NAO_MO.data, c='black',linewidth=1,  label='MONSOON')
 plt.fill_between(times, 0, np.ma.masked_where(NAO_MO.data <= 0, NAO_MO.data) , alpha=0.5, facecolor='red')
 plt.fill_between(times, 0, np.ma.masked_where(NAO_MO.data >= 0, NAO_MO.data) , alpha=0.5, facecolor='blue')
 plt.xlabel('Time(years)', fontsize=14)
 plt.ylabel('NAO Index', fontsize=14)
 plt.title(title_2, fontsize=16)
 #plt.ticklabel_format(axis='y',useOffset=False) 
 #ax.xaxis.set_major_locator(mdates.YearLocator(20))
 plt.legend()
 
 #iplt.show()
 return 

##########################################################################################################################################



##loading  monthly data 
AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_SLP_18502050_monthly.nc')
MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_SLP_18502050_monthly.nc')

NAO_AR = NAO(AR, seasonal = 1)
NAO_MO = NAO(MO, seasonal = 1)
Plot(NAO_MO, NAO_AR, title_1= 'Winter (DJFM) NAO Index ARCHER 1850-2050', title_2 = 'Winter (DJFM) NAO Index MONSOON 1850-2050' )


##loading annual data 
AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_SLP_U_V_18502050.nc',  iris.Constraint('air_pressure_at_sea_level'))
MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_SLP_U_V_18502050.nc',  iris.Constraint('air_pressure_at_sea_level'))

NAO_AR = NAO(AR, seasonal = 0)
NAO_MO = NAO(MO, seasonal = 0)
Plot(NAO_MO, NAO_AR, title_1= 'Annual NAO Index ARCHER 1850-2050', title_2 = 'Annual NAO Index MONSOON 1850-2050' )

iplt.show()
