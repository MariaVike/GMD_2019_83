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

#from shiftedColorMap import*
from scipy import stats

'''
############################################## CALCULATE ENSO INDEX #############################################################

#MO_SST_cubes=iris.load('/gws/nopw/j04/pmip4_vol1/users/vittoria/seaice_monthly_uar766_*.nc', iris.Constraint('sea surface temperature', latitude= lambda lat: -5 <= lat <= 5, longitude= lambda lon: -170 <= lon <= -120   ))

LIG_SST_cubes=iris.load('/gws/nopw/j04/pmip4_vol1/users/vittoria/seaice_monthly_uba937_*.nc',  iris.Constraint('sea surface temperature', latitude= lambda lat: -5 <= lat <= 5, longitude= lambda lon: -170 <= lon <= -120  ))

#concatenate cubes in one single dataset to plot
#iris.util.unify_time_units(MO_SST_cubes)
iris.util.unify_time_units(LIG_SST_cubes)

#MO_SST=MO_SST_cubes.concatenate_cube()
LIG_SST=LIG_SST_cubes.concatenate_cube()

#qplt.contourf(LIG_SST[0,:,:])
#plt.show()

LIG_SST.coord('latitude').guess_bounds()
LIG_SST.coord('longitude').guess_bounds()
area_weights_LIG= iris.analysis.cartography.area_weights(LIG_SST)

#MO_SST.coord('latitude').guess_bounds()
#MO_SST.coord('longitude').guess_bounds()
#area_weights_MO= iris.analysis.cartography.area_weights(MO_SST)


LIG_SST_av = LIG_SST.collapsed(['latitude', 'longitude'],
                                 iris.analysis.MEAN,
                                 weights=area_weights_LIG)

#MO_SST_av = MO_SST.collapsed(['latitude', 'longitude'],
                                 #iris.analysis.MEAN,
                                # weights=area_weights_MO)

#print LIG_SST_av
#print LIG_SST_av.coords('latitude')
#print LIG_SST_av.coords('longitude')

#calculate the mean over a 30-year reference period: 
LIG_SST_tmean=LIG_SST_av[:30*12].collapsed('time', iris.analysis.MEAN)
#MO_SST_tmean=MO_SST_av[:30*12].collapsed('time', iris.analysis.MEAN)

#print LIG_SST
#print MO_SST

#calculate anomalies:
ENSO_LIG= LIG_SST_av-LIG_SST_tmean
#ENSO_MO= MO_SST_av-MO_SST_tmean


#filter the time series with a 3-month running mean:
months_LIG=ENSO_LIG.shape[0]
ENSO_LIG_3month=iris.cube.CubeList() 

for i in range (0, months_LIG-3):
   ENSO_LIG_3month.append(ENSO_LIG[i:i+3].collapsed('time', iris.analysis.MEAN))
  
ENSO_LIG_3month= ENSO_LIG_3month.merge_cube()
print ENSO_LIG_3month

#months_MO=ENSO_MO.shape[0]
#ENSO_MO_3month=iris.cube.CubeList() 

#for i in range (0, months_MO-3):
 #   ENSO_MO_3month.append(ENSO_MO[i:i+3].collapsed('time', iris.analysis.MEAN)) 

#ENSO_MO_3month= ENSO_MO_3month.merge_cube()
#print ENSO_MO_3month


iris.save(ENSO_LIG_3month, '/gws/nopw/j04/pmip4_vol1/users/vittoria/ENSO_LIG_18502050.nc')
#iris.save(ENSO_MO_3month, 'ENSO_MO_18502050.nc')
'''

############################################## PLOTTING #############################################################

ENSO_LIG_3month=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/ENSO_LIG_18502050.nc')
ENSO_MO_3month=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/ENSO_MO_18502050.nc')

##select interval
#ENSO_LIG_3month=ENSO_LIG_3month[:50*12]
#ENSO_MO_3month=ENSO_MO_3month[:50*12]

mean_LIG=ENSO_LIG_3month.collapsed('time', iris.analysis.MEAN)
stdev_LIG=ENSO_LIG_3month.collapsed('time', iris.analysis.STD_DEV)
mean_MO=ENSO_MO_3month.collapsed('time', iris.analysis.MEAN)
stdev_MO=ENSO_MO_3month.collapsed('time', iris.analysis.STD_DEV)
print mean_LIG.data, mean_MO.data, stdev_LIG.data, stdev_MO.data

months=ENSO_LIG_3month.shape[0]
times= pd.date_range(start='1850-01-01', periods=months, freq='MS') #create array of dates with pandas 
##edit times for select interval
#times=times[:50*12]

nino_thrs=np.empty(ENSO_LIG_3month.shape[0]); nino_thrs.fill(0.5)
nina_thrs=np.empty(ENSO_LIG_3month.shape[0]); nina_thrs.fill(-0.5)

fig, ax = plt.subplots(1)
fig.autofmt_xdate()
plt.plot(times, ENSO_LIG_3month.data, c='black', linestyle= '--',linewidth=2, label='LIG')
plt.plot(times, ENSO_MO_3month.data, c='grey',linewidth=2,  label='PI')
plt.plot(times, nino_thrs, c='red' )
plt.plot(times, nina_thrs, c='blue')
plt.xlabel('Date', fontsize=14)
plt.ylabel('3.4 ENSO Index', fontsize=14)
plt.title(r'3.4 ENSO Index LIG vs PI 1850-2050', fontsize=16)
plt.ticklabel_format(axis='y',useOffset=False)  #only on y axis if we are using a date array created with pandas
ax.xaxis.set_major_locator(mdates.YearLocator(10))
plt.legend()


#scatter plot
fig, ax = plt.subplots(1)
Nino_LIG=np.ma.masked_where(ENSO_LIG_3month.data < 0.5, ENSO_LIG_3month.data)
Nina_LIG=np.ma.masked_where(ENSO_LIG_3month.data > -0.5, ENSO_LIG_3month.data)
Nino_MO=np.ma.masked_where(ENSO_MO_3month.data < 0.5, ENSO_MO_3month.data)
Nina_MO=np.ma.masked_where(ENSO_MO_3month.data > -0.5, ENSO_MO_3month.data)

print Nino_LIG.count(), Nino_MO.count() #print the number of unmasked values 
print Nina_LIG.count(), Nina_MO.count()

line=np.arange(-4,4,1)
plt.plot(line, line, c='black')
plt.scatter(ENSO_MO_3month.data, ENSO_LIG_3month.data, c='grey')
agree=plt.scatter(Nino_MO, Nino_LIG, c='red') #MO & AR agree : NINO for both 
plt.scatter(Nina_MO, Nina_LIG, c='red') #MO & AR agree : NINA for both 
disagree=plt.scatter(Nino_MO, Nina_LIG, c='blue') # MO & AR disagree
plt.scatter(Nina_MO, Nino_LIG, c='blue') # MO & AR disagree
plt.xlabel('ENSO Index PI', fontsize=14)
plt.ylabel('ENSO Index LIG', fontsize=14)
plt.title(r'Scatter Plot - ENSO Index LIG vs PI 1850-2050', fontsize=16)
plt.legend((agree, disagree), ('PI & LIG agree', 'PI & LIG disagree'), loc=1, fontsize=13)


####second version of same plot, LIG and PI on 2 seprate sub-plots #####################

fig = plt.figure()

ax = fig.add_subplot(211)

fig.autofmt_xdate()
plt.plot(times, ENSO_LIG_3month.data, c='black',linewidth=1, label='LIG')
plt.fill_between(times, 0.5, np.ma.masked_where(ENSO_LIG_3month.data <= 0.5, ENSO_LIG_3month.data) , alpha=0.5, facecolor='orangered') 
plt.fill_between(times, -0.5, np.ma.masked_where(ENSO_LIG_3month.data >= -0.5, ENSO_LIG_3month.data) , alpha=0.5, facecolor='darkgreen')
plt.xlabel('Time(years)', fontsize=14)
plt.ylabel('ENSO Index', fontsize=14)
plt.title('ENSO Index LIG 1850-2050', fontsize=16)
plt.ticklabel_format(axis='y',useOffset=False) 
ax.xaxis.set_major_locator(mdates.YearLocator(10))
plt.legend()


ax = fig.add_subplot(212)

fig.autofmt_xdate()
plt.plot(times, ENSO_MO_3month.data, c='black',linewidth=1,  label='MONSOON')
plt.fill_between(times, 0.5, np.ma.masked_where(ENSO_MO_3month.data <= 0.5, ENSO_MO_3month.data) , alpha=0.5, facecolor='orangered')
plt.fill_between(times, -0.5, np.ma.masked_where(ENSO_MO_3month.data >= -0.5, ENSO_MO_3month.data) , alpha=0.5, facecolor='darkgreen')
plt.xlabel('Time(years)', fontsize=14)
plt.ylabel('ENSO Index', fontsize=14)
plt.title('ENSO Index PI 1850-2050', fontsize=16)
plt.ticklabel_format(axis='y',useOffset=False) 
ax.xaxis.set_major_locator(mdates.YearLocator(10))
plt.legend()

iplt.show()



