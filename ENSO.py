import iris 
import iris.plot as iplt
import iris.quickplot as qplt

import iris.analysis.cartography 
import cartopy.crs as ccrs
import cartopy.feature as cfe
from mpl_toolkits.basemap import Basemap, maskoceans

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from iris.experimental.equalise_cubes import *
from iris.experimental.animate import*

import warnings
import glob
import os 

import pandas as pd
import scipy
import matplotlib.dates as mdates

from shiftedColorMap import*
from scipy import stats

'''
############################################## CALCULATE ENSO INDEX #############################################################

MO_SST_cubes=iris.load('/nerc/n02/n02/vittoria/seaice_monthly_uar766_*.nc', iris.Constraint('sea surface temperature', latitude= lambda lat: -5 <= lat <= 5, longitude= lambda lon: -170 <= lon <= -120   ))

AR_SST_cubes=iris.load('/nerc/n02/n02/vittoria/seaice_monthly_uas245_*.nc',  iris.Constraint('sea surface temperature', latitude= lambda lat: -5 <= lat <= 5, longitude= lambda lon: -170 <= lon <= -120  ))

#concatenate cubes in one single dataset to plot
iris.util.unify_time_units(MO_SST_cubes)
iris.util.unify_time_units(AR_SST_cubes)

MO_SST=MO_SST_cubes.concatenate_cube()
AR_SST=AR_SST_cubes.concatenate_cube()

#qplt.contourf(AR_SST[0,:,:])
#plt.show()

AR_SST.coord('latitude').guess_bounds()
AR_SST.coord('longitude').guess_bounds()
area_weights_AR= iris.analysis.cartography.area_weights(AR_SST)

MO_SST.coord('latitude').guess_bounds()
MO_SST.coord('longitude').guess_bounds()
area_weights_MO= iris.analysis.cartography.area_weights(MO_SST)


AR_SST_av = AR_SST.collapsed(['latitude', 'longitude'],
                                 iris.analysis.MEAN,
                                 weights=area_weights_AR)

MO_SST_av = MO_SST.collapsed(['latitude', 'longitude'],
                                 iris.analysis.MEAN,
                                 weights=area_weights_MO)

#print AR_SST_av
#print AR_SST_av.coords('latitude')
#print AR_SST_av.coords('longitude')

#calculate the mean over a 30-year reference period: 
AR_SST_tmean=AR_SST_av[:30*12].collapsed('time', iris.analysis.MEAN)
MO_SST_tmean=MO_SST_av[:30*12].collapsed('time', iris.analysis.MEAN)

#print AR_SST
#print MO_SST

#calculate anomalies:
ENSO_AR= AR_SST_av-AR_SST_tmean
ENSO_MO= MO_SST_av-MO_SST_tmean


#filter the time series with a 3-month running mean:
months_AR=ENSO_AR.shape[0]
ENSO_AR_3month=iris.cube.CubeList() 

for i in range (0, months_AR-3):
   ENSO_AR_3month.append(ENSO_AR[i:i+3].collapsed('time', iris.analysis.MEAN))
  
ENSO_AR_3month= ENSO_AR_3month.merge_cube()
print ENSO_AR_3month


months_MO=ENSO_MO.shape[0]
ENSO_MO_3month=iris.cube.CubeList() 

for i in range (0, months_MO-3):
    ENSO_MO_3month.append(ENSO_MO[i:i+3].collapsed('time', iris.analysis.MEAN)) 

ENSO_MO_3month= ENSO_MO_3month.merge_cube()
print ENSO_MO_3month


iris.save(ENSO_AR_3month, 'ENSO_AR_18502050.nc')
iris.save(ENSO_MO_3month, 'ENSO_MO_18502050.nc')
'''

############################################## PLOTTING #############################################################

ENSO_AR_3month=iris.load_cube('/nerc/n02/n02/vittoria/ENSO_AR_18502050.nc')
ENSO_MO_3month=iris.load_cube('/nerc/n02/n02/vittoria/ENSO_MO_18502050.nc')

##select interval
#ENSO_AR_3month=ENSO_AR_3month[:50*12]
#ENSO_MO_3month=ENSO_MO_3month[:50*12]

mean_AR=ENSO_AR_3month.collapsed('time', iris.analysis.MEAN)
stdev_AR=ENSO_AR_3month.collapsed('time', iris.analysis.STD_DEV)
mean_MO=ENSO_MO_3month.collapsed('time', iris.analysis.MEAN)
stdev_MO=ENSO_MO_3month.collapsed('time', iris.analysis.STD_DEV)
print mean_AR.data, mean_MO.data, stdev_AR.data, stdev_MO.data

months=ENSO_AR_3month.shape[0]
times= pd.date_range(start='1850-01-01', periods=months, freq='MS') #create array of dates with pandas 
##edit times for select interval
#times=times[:50*12]

nino_thrs=np.empty(ENSO_AR_3month.shape[0]); nino_thrs.fill(0.5)
nina_thrs=np.empty(ENSO_AR_3month.shape[0]); nina_thrs.fill(-0.5)

fig, ax = plt.subplots(1)
fig.autofmt_xdate()
plt.plot(times, ENSO_AR_3month.data, c='black', linestyle= '--',linewidth=2, label='ARCHER')
plt.plot(times, ENSO_MO_3month.data, c='grey',linewidth=2,  label='MONSOON')
plt.plot(times, nino_thrs, c='red' )
plt.plot(times, nina_thrs, c='blue')
plt.xlabel('Date', fontsize=14)
plt.ylabel('3.4 ENSO Index', fontsize=14)
plt.title(r'3.4 ENSO Index AR vs MO 1850-2050', fontsize=16)
plt.ticklabel_format(axis='y',useOffset=False)  #only on y axis if we are using a date array created with pandas
ax.xaxis.set_major_locator(mdates.YearLocator(10))
plt.legend()


#scatter plot
fig, ax = plt.subplots(1)
Nino_AR=np.ma.masked_where(ENSO_AR_3month.data < 0.5, ENSO_AR_3month.data)
Nina_AR=np.ma.masked_where(ENSO_AR_3month.data > -0.5, ENSO_AR_3month.data)
Nino_MO=np.ma.masked_where(ENSO_MO_3month.data < 0.5, ENSO_MO_3month.data)
Nina_MO=np.ma.masked_where(ENSO_MO_3month.data > -0.5, ENSO_MO_3month.data)

print Nino_AR.count(), Nino_MO.count() #print the number of unmasked values 
print Nina_AR.count(), Nina_MO.count()

line=np.arange(-4,4,1)
plt.plot(line, line, c='black')
plt.scatter(ENSO_MO_3month.data, ENSO_AR_3month.data, c='grey')
agree=plt.scatter(Nino_MO, Nino_AR, c='red') #MO & AR agree : NINO for both 
plt.scatter(Nina_MO, Nina_AR, c='red') #MO & AR agree : NINA for both 
disagree=plt.scatter(Nino_MO, Nina_AR, c='blue') # MO & AR disagree
plt.scatter(Nina_MO, Nino_AR, c='blue') # MO & AR disagree
plt.xlabel('ENSO Index MO', fontsize=14)
plt.ylabel('ENSO Index AR', fontsize=14)
plt.title(r'Scatter Plot - ENSO Index AR vs MO 1850-2050', fontsize=16)
plt.legend((agree, disagree), ('MO & AR agree', 'MO & AR disagree'), loc=1, fontsize=13)


####second version of same plot, AR and MO on 2 seprate sub-plots #####################

fig = plt.figure()

ax = fig.add_subplot(211)

fig.autofmt_xdate()
plt.plot(times, ENSO_AR_3month.data, c='black',linewidth=1, label='ARCHER')
plt.fill_between(times, 0.5, np.ma.masked_where(ENSO_AR_3month.data <= 0.5, ENSO_AR_3month.data) , alpha=0.5, facecolor='orangered') 
plt.fill_between(times, -0.5, np.ma.masked_where(ENSO_AR_3month.data >= -0.5, ENSO_AR_3month.data) , alpha=0.5, facecolor='darkgreen')
plt.xlabel('Time(years)', fontsize=14)
plt.ylabel('SAM Index', fontsize=14)
plt.title('ENSO Index ARCHER 1850-2050', fontsize=16)
plt.ticklabel_format(axis='y',useOffset=False) 
ax.xaxis.set_major_locator(mdates.YearLocator(10))
plt.legend()


ax = fig.add_subplot(212)

fig.autofmt_xdate()
plt.plot(times, ENSO_MO_3month.data, c='black',linewidth=1,  label='MONSOON')
plt.fill_between(times, 0.5, np.ma.masked_where(ENSO_MO_3month.data <= 0.5, ENSO_MO_3month.data) , alpha=0.5, facecolor='orangered')
plt.fill_between(times, -0.5, np.ma.masked_where(ENSO_MO_3month.data >= -0.5, ENSO_MO_3month.data) , alpha=0.5, facecolor='darkgreen')
plt.xlabel('Time(years)', fontsize=14)
plt.ylabel('SAM Index', fontsize=14)
plt.title('ENSO Index MONSOON 1850-2050', fontsize=16)
plt.ticklabel_format(axis='y',useOffset=False) 
ax.xaxis.set_major_locator(mdates.YearLocator(10))
plt.legend()

iplt.show()



