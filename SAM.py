import iris 
import iris.plot as iplt
import iris.quickplot as qplt
import iris.coord_categorisation
from iris.util import *

import iris.analysis.cartography 
import cartopy.crs as ccrs
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


#################################### SAM INDEX ARCHER ######################################################################

#P40_AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_SLP_18502050_monthly.nc', iris.Constraint(latitude=lambda l: l == -40.625 ))
#P65_AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_SLP_18502050_monthly.nc', iris.Constraint(latitude=lambda l: l == -65.625 ))      
P40_AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_SLP_U_V_18502050.nc',  iris.Constraint('air_pressure_at_sea_level', latitude=lambda l: l == -40.625 ))
P65_AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_SLP_U_V_18502050.nc', iris.Constraint('air_pressure_at_sea_level', latitude=lambda l: l == -65.625 ))    

#calcualate  zonal mean
P40_znm_AR=P40_AR.collapsed('longitude', iris.analysis.MEAN)
P65_znm_AR=P65_AR.collapsed('longitude', iris.analysis.MEAN)

#mean over first 30 years to be used to calculate anomalies
### monthly data ######
#P40_mean_AR=P40_znm_AR[:30*12].collapsed('time', iris.analysis.MEAN)
#P40_stdev_AR=P40_znm_AR[:30*12].collapsed('time', iris.analysis.STD_DEV)

#P65_mean_AR=P65_znm_AR[:30*12].collapsed('time', iris.analysis.MEAN)
#P65_stdev_AR=P65_znm_AR[:30*12].collapsed('time', iris.analysis.STD_DEV)

### annual data ######
P40_mean_AR=P40_znm_AR[:30].collapsed('time', iris.analysis.MEAN)
P40_stdev_AR=P40_znm_AR[:30].collapsed('time', iris.analysis.STD_DEV)

P65_mean_AR=P65_znm_AR[:30].collapsed('time', iris.analysis.MEAN)
P65_stdev_AR=P65_znm_AR[:30].collapsed('time', iris.analysis.STD_DEV)

#normalize each month separately
P40S_AR=(P40_znm_AR-P40_mean_AR)/P40_stdev_AR
P65S_AR=(P65_znm_AR-P65_mean_AR)/P65_stdev_AR

SAM_AR=P40S_AR - P65S_AR

#################################### SAM INDEX Met Office #########################################################################################################

#P40_MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_SLP_18502050_monthly.nc', iris.Constraint(latitude=lambda l: l == -40.625 ))
#P65_MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_SLP_18502050_monthly.nc', iris.Constraint(latitude=lambda l: l == -65.625 ))
P40_MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_SLP_U_V_18502050.nc', iris.Constraint('air_pressure_at_sea_level', latitude=lambda l: l == -40.625 ))
P65_MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_SLP_U_V_18502050.nc', iris.Constraint('air_pressure_at_sea_level', latitude=lambda l: l == -65.625 ))

#calcualate  zonal mean
P40_znm_MO=P40_MO.collapsed('longitude', iris.analysis.MEAN)
P65_znm_MO=P65_MO.collapsed('longitude', iris.analysis.MEAN)

#mean over first 30 years to be used to calculate anomalies
P40_mean_MO=P40_znm_MO[:30*12].collapsed('time', iris.analysis.MEAN)
P40_stdev_MO=P40_znm_MO[:30*12].collapsed('time', iris.analysis.STD_DEV)

P65_mean_MO=P65_znm_MO[:30*12].collapsed('time', iris.analysis.MEAN)
P65_stdev_MO=P65_znm_MO[:30*12].collapsed('time', iris.analysis.STD_DEV)


P40S_MO=(P40_znm_MO-P40_mean_MO)/P40_stdev_MO
P65S_MO=(P65_znm_MO-P65_mean_MO)/P65_stdev_MO

SAM_MO=P40S_MO - P65S_MO

'''
########################### 3-month running mean ###################################################################################################################

months_AR=SAM_AR.shape[0]
SAM_AR_3month=iris.cube.CubeList() 

for i in range (0, months_AR-3):
   SAM_AR_3month.append(SAM_AR[i:i+3].collapsed('time', iris.analysis.MEAN))
  
SAM_AR_3month= SAM_AR_3month.merge_cube()
print SAM_AR_3month


months_MO=SAM_MO.shape[0]
SAM_MO_3month=iris.cube.CubeList() 

for i in range (0, months_MO-3):
    SAM_MO_3month.append(SAM_MO[i:i+3].collapsed('time', iris.analysis.MEAN)) 

SAM_MO_3month= SAM_MO_3month.merge_cube()
print SAM_MO_3month


iris.save(SAM_AR_3month, 'SAM_AR_18502050_3month_runn.nc')
iris.save(SAM_MO_3month, 'SAM_MO_18502050_3month_runn.nc')
'''

################################################## PLOt SAM INDEXES ###########################################################################################

#SAM_AR=iris.load_cube('/nerc/n02/n02/vittoria/SAM_AR_18502050_3month_runn.nc')
#SAM_MO=iris.load_cube('/nerc/n02/n02/vittoria/SAM_MO_18502050_3month_runn.nc')

##select interval
#SAM_AR=SAM_AR[:50*12]
#SAM_MO=SAM_MO[:50*12]

mean_AR=SAM_AR.collapsed('time', iris.analysis.MEAN)
stdev_AR=SAM_AR.collapsed('time', iris.analysis.STD_DEV)
mean_MO=SAM_MO.collapsed('time', iris.analysis.MEAN)
stdev_MO=SAM_MO.collapsed('time', iris.analysis.STD_DEV)
print mean_AR.data, mean_MO.data, stdev_AR.data, stdev_MO.data

months=SAM_AR.shape[0]
times= pd.date_range(start='1850-01-01', periods=months, freq='AS') #for annual values
#times= pd.date_range(start='1850-01-01', periods=months, freq='MS') #for monthly values
##edit times for select interval
#times=times[:50*12]
times=np.arange(0,months,1) #only number of show years from 1850 9for annual SAM)

zero_line=np.empty(SAM_AR.shape[0]); zero_line.fill(0)

fig = plt.figure()

ax = fig.add_subplot(211)

fig.autofmt_xdate()
plt.plot(times, zero_line, c='black' )
plt.plot(times, SAM_AR.data, c='black',linewidth=1, label='ARCHER')
plt.fill_between(times, 0, np.ma.masked_where(SAM_AR.data <= 0, SAM_AR.data) , alpha=0.5, facecolor='orange')
plt.fill_between(times, 0, np.ma.masked_where(SAM_AR.data >= 0, SAM_AR.data) , alpha=0.5, facecolor='blue')
plt.xlabel('Time(years)', fontsize=14)
plt.ylabel('SAM Index', fontsize=14)
plt.title(r'SAM Index AR 1850-2050', fontsize=16)
#plt.ticklabel_format(useOffset=False)
#plt.ticklabel_format(axis='y',useOffset=False)  #only on y axis if we are using a date array created with pandas
#ax.xaxis.set_major_locator(mdates.YearLocator(20))
plt.legend()


ax = fig.add_subplot(212)

fig.autofmt_xdate()
plt.plot(times, zero_line, c='black' )
plt.plot(times, SAM_MO.data, c='black',linewidth=1,  label='MONSOON')
plt.fill_between(times, 0, np.ma.masked_where(SAM_MO.data <= 0, SAM_MO.data) , alpha=0.5, facecolor='orange')
plt.fill_between(times, 0, np.ma.masked_where(SAM_MO.data >= 0, SAM_MO.data) , alpha=0.5, facecolor='blue')
plt.xlabel('Time(years)', fontsize=14)
plt.ylabel('SAM Index', fontsize=14)
plt.title(r'SAM Index MO 1850-2050', fontsize=16)
#plt.ticklabel_format(useOffset=False)
#plt.ticklabel_format(axis='y',useOffset=False)  #only on y axis if we are using a date array created with pandas
#ax.xaxis.set_major_locator(mdates.YearLocator(20))
plt.legend()

iplt.show()



###################### PLOT U at 850 hPa for SH JET POSITION AND STRENGTH ######################################

U850_AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_U_V_3D_18502050.nc', iris.Constraint('x_wind', pressure=lambda p: p == 850))      
U850_MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_U_V_3D_18502050.nc', iris.Constraint('x_wind', pressure=lambda p: p == 850))
#computing 200yr means for AR and MO
U850_MO_mean=U850_MO.collapsed('time', iris.analysis.MEAN)
U850_AR_mean=U850_AR.collapsed('time', iris.analysis.MEAN)
#computing mean of MO-AR diff and STdev of MO for S/N
U850_mean_diff= (U850_MO[:,:,:]-U850_AR[:,:,:]).collapsed('time', iris.analysis.MEAN)
U850_MO_stdev=U850_MO.collapsed('time', iris.analysis.STD_DEV)


fig=plt.figure(figsize=(20,5), tight_layout=True)
fig.suptitle('U at 850 hPa 200-year mean (m/s)' , x=0.5, y=0.05, fontsize=28)

ax = fig.add_subplot(141)
ax=plt.gca()
bmap=Basemap(projection='ortho',lat_0=-90,lon_0=0,resolution='l', ax=ax)
#map.shadedrelief(scale=0.25)

#get the lon,lat values from cube in order to plot: 
#lon=U850_MO_mean.coord('longitude').points
lon=np.linspace(0, 360, 192) 
lat=U850_MO_mean.coord('latitude').points
x,y=bmap(*np.meshgrid(lon,lat))

#orig_cmap=matplotlib.cm.RdBu_r 
#shifted_cmap=shiftedColorMap(orig_cmap, midpoint=0.42, name='shifted')
levels=40
contours=bmap.contourf(x,y, U850_MO_mean[:,:].data, levels, cmap='RdBu_r') 

bmap.drawcoastlines()
bmap.drawparallels(np.arange(-80.,81.,20.))
bmap.drawmeridians(np.arange(-180.,181.,20.))
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label('U at 850 hPa', fontsize=16) 
plt.title(r'MONSOON', fontsize=16)

ax = fig.add_subplot(142)
ax=plt.gca()
bmap=Basemap(projection='ortho',lat_0=-90,lon_0=0,resolution='l', ax=ax)
#map.shadedrelief(scale=0.25)

#get the lon,lat values from cube in order to plot: 
lon=np.linspace(0, 360, 192) 
lat=U850_AR_mean.coord('latitude').points
x,y=bmap(*np.meshgrid(lon,lat))

#orig_cmap=matplotlib.cm.RdBu_r 
#shifted_cmap=shiftedColorMap(orig_cmap, midpoint=0.42, name='shifted')
levels=40
contours=bmap.contourf(x,y, U850_AR_mean[:,:].data, levels, cmap='RdBu_r') 

bmap.drawcoastlines()
bmap.drawparallels(np.arange(-80.,81.,20.))
bmap.drawmeridians(np.arange(-180.,181.,20.))
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label('U at 850 hPa', fontsize=16)
plt.title(r'ARCHER', fontsize=16)


ax = fig.add_subplot(143)
ax=plt.gca()
bmap=Basemap(projection='ortho',lat_0=-90,lon_0=0,resolution='l', ax=ax)
#map.shadedrelief(scale=0.25)

x,y=bmap(*np.meshgrid(lon,lat))
contours=bmap.contourf(x,y, (U850_mean_diff.data ) , levels, cmap='PuOr_r') 

bmap.drawcoastlines()
bmap.drawparallels(np.arange(-80.,81.,20.))
bmap.drawmeridians(np.arange(-180.,181.,20.))
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label('U$_{MO}$ - U$_{AR}$ at 850 hPa', fontsize=16) 
plt.title(r'MONSOON - ARCHER', fontsize=16)


ax = fig.add_subplot(144)
ax=plt.gca()
bmap=Basemap(projection='ortho',lat_0=-90,lon_0=0,resolution='l', ax=ax)
#map.shadedrelief(scale=0.25)
x,y=bmap(*np.meshgrid(lon,lat))

U850_MO_stdev.data=np.ma.masked_less_equal(U850_MO_stdev.data, 0.01)
S_to_N=U850_mean_diff/U850_MO_stdev
contours=bmap.contourf(x,y, abs(S_to_N.data) , levels, cmap='GnBu_r') 
#contours=bmap.contourf(x,y, U850_MO_stdev.data , levels, cmap='GnBu_r') 

bmap.drawcoastlines()
bmap.drawparallels(np.arange(-80.,81.,20.))
bmap.drawmeridians(np.arange(-180.,181.,20.))
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label('S/N of U$_{MO}$ - U$_{AR}$ at 850 hPa', fontsize=16) 
plt.title(r'S/N of MONSOON - ARCHER', fontsize=16)

#plt.show()

###################### PLOT zonal mean cross section of U vs height/pressure  ######################################

U_AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_U_V_3D_18502050.nc', iris.Constraint('x_wind', latitude= lambda lat: -90 <= lat <= 0))       
U_MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_U_V_3D_18502050.nc', iris.Constraint('x_wind', latitude= lambda lat: -90 <= lat <= 0))

U_MO_zmn=U_MO.collapsed('longitude', iris.analysis.MEAN)
U_AR_zmn=U_AR.collapsed('longitude', iris.analysis.MEAN)

U_MO_mean=U_MO_zmn.collapsed('time', iris.analysis.MEAN)
U_AR_mean=U_AR_zmn.collapsed('time', iris.analysis.MEAN)
#print U_MO_mean

#computing mean of MO-AR diff and STdev of MO for S/N
U_mean_diff= (U_MO_zmn[:,:,:]-U_AR_zmn[:,:,:]).collapsed('time', iris.analysis.MEAN)
U_MO_stdev=U_MO_zmn.collapsed('time', iris.analysis.STD_DEV)


fig=plt.figure(figsize=(20,5), tight_layout=True)
fig.suptitle(r'Zonal Mean Zonal Wind $\bar{U}$ 200-year mean (m/s)' , x=0.5, y=0.05, fontsize=28)

ax = fig.add_subplot(141)
ax=plt.gca()
contours=iplt.contourf(U_MO_mean[:,:],cmap='RdBu_r') #,levels)
#contours_2=iplt.contour(U_MO_mean[:,:], levels, colors='black', linewidths=1)
#plt.clabel(contours_2)
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label(r'$\bar{U}_{MO}$', fontsize=16) 
plt.title(r'MONSOON', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
plt.ylabel('Height (hPa)', fontsize=14)
plt.xticks([-80, -60, -40, -20, 0])


ax = fig.add_subplot(142)
ax=plt.gca()
contours=iplt.contourf(U_AR_mean[:,:],cmap='RdBu_r') #,levels
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label(r'$\bar{U}_{AR}$', fontsize=16) 
plt.title(r'ARCHER', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
#plt.ylabel('Height (hPa)', fontsize=14)
plt.xticks([-80, -60, -40, -20, 0])

ax = fig.add_subplot(143)
ax=plt.gca()
contours=iplt.contourf(U_mean_diff ,cmap='PuOr_r') #,levels
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label(r'$\bar{U}_{MO} - \bar{U}_{AR}$', fontsize=16) 
plt.title(r'MONSOON - ARCHER ', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
#plt.ylabel('Height (hPa)', fontsize=14)
plt.xticks([-80, -60, -40, -20, 0])



ax = fig.add_subplot(144)
ax=plt.gca()

S_to_N=(U_mean_diff/U_MO_stdev)
S_to_N.data=abs(S_to_N.data)
contours=iplt.contourf(S_to_N ,cmap='GnBu_r')
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label(r'$S/N of \bar{U}_{MO} - \bar{U}_{AR}$', fontsize=16) 
plt.title(r'S/N of MONSOON - ARCHER ', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
#plt.ylabel('Height (hPa)', fontsize=14)
plt.xticks([-80, -60, -40, -20, 0])


plt.show()


