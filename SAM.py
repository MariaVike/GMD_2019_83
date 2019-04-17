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


#################################### SAM INDEX PI ######################################################################

#P40_PI=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/uar766_SLP_18502050_monthly.nc', iris.Constraint(latitude=lambda l: l == -40.625 ))
#P65_PI=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/uar766_SLP_18502050_monthly.nc', iris.Constraint(latitude=lambda l: l == -65.625 ))      
P40_PI=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/uar766_SLP_U_V_18502050.nc',  iris.Constraint('air_pressure_at_sea_level', latitude=lambda l: l == -40.625 ))
P65_PI=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/uar766_SLP_U_V_18502050.nc', iris.Constraint('air_pressure_at_sea_level', latitude=lambda l: l == -65.625 ))    

#calcualate  zonal mean
P40_znm_PI=P40_PI.collapsed('longitude', iris.analysis.MEAN)
P65_znm_PI=P65_PI.collapsed('longitude', iris.analysis.MEAN)

#mean over first 30 years to be used to calculate anomalies
### monthly data ######
#P40_mean_PI=P40_znm_PI[:30*12].collapsed('time', iris.analysis.MEAN)
#P40_stdev_PI=P40_znm_PI[:30*12].collapsed('time', iris.analysis.STD_DEV)

#P65_mean_PI=P65_znm_PI[:30*12].collapsed('time', iris.analysis.MEAN)
#P65_stdev_PI=P65_znm_PI[:30*12].collapsed('time', iris.analysis.STD_DEV)

### annual data ######
P40_mean_PI=P40_znm_PI[:30].collapsed('time', iris.analysis.MEAN)
P40_stdev_PI=P40_znm_PI[:30].collapsed('time', iris.analysis.STD_DEV)

P65_mean_PI=P65_znm_PI[:30].collapsed('time', iris.analysis.MEAN)
P65_stdev_PI=P65_znm_PI[:30].collapsed('time', iris.analysis.STD_DEV)

#normalize each month separately
P40S_PI=(P40_znm_PI-P40_mean_PI)/P40_stdev_PI
P65S_PI=(P65_znm_PI-P65_mean_PI)/P65_stdev_PI

SAM_PI=P40S_PI - P65S_PI

#################################### SAM INDEX LIG #########################################################################################################

#P40_LIG=iris.load_cube('/nerc/n02/n02/vittoria/uar766_SLP_18502050_monthly.nc', iris.Constraint(latitude=lambda l: l == -40.625 ))
#P65_LIG=iris.load_cube('/nerc/n02/n02/vittoria/uar766_SLP_18502050_monthly.nc', iris.Constraint(latitude=lambda l: l == -65.625 ))
P40_LIG=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/uba937_SLP_U_V_18502050.nc', iris.Constraint('air_pressure_at_sea_level', latitude=lambda l: l == -40.625 ))
P65_LIG=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/uba937_SLP_U_V_18502050.nc', iris.Constraint('air_pressure_at_sea_level', latitude=lambda l: l == -65.625 ))

#calcualate  zonal mean
P40_znm_LIG=P40_LIG.collapsed('longitude', iris.analysis.MEAN)
P65_znm_LIG=P65_LIG.collapsed('longitude', iris.analysis.MEAN)

#mean over first 30 years to be used to calculate anomalies
P40_mean_LIG=P40_znm_LIG[:30*12].collapsed('time', iris.analysis.MEAN)
P40_stdev_LIG=P40_znm_LIG[:30*12].collapsed('time', iris.analysis.STD_DEV)

P65_mean_LIG=P65_znm_LIG[:30*12].collapsed('time', iris.analysis.MEAN)
P65_stdev_LIG=P65_znm_LIG[:30*12].collapsed('time', iris.analysis.STD_DEV)


P40S_LIG=(P40_znm_LIG-P40_mean_LIG)/P40_stdev_LIG
P65S_LIG=(P65_znm_LIG-P65_mean_LIG)/P65_stdev_LIG

SAM_LIG=P40S_LIG - P65S_LIG

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

#SAM_PI=iris.load_cube('/nerc/n02/n02/vittoria/SAM_PI_18502050_3month_runn.nc')
#SAM_LIG=iris.load_cube('/nerc/n02/n02/vittoria/SAM_LIG_18502050_3month_runn.nc')

##select interval
#SAM_PI=SAM_PI[:50*12]
#SAM_LIG=SAM_LIG[:50*12]

mean_PI=SAM_PI.collapsed('time', iris.analysis.MEAN)
stdev_PI=SAM_PI.collapsed('time', iris.analysis.STD_DEV)
mean_LIG=SAM_LIG.collapsed('time', iris.analysis.MEAN)
stdev_LIG=SAM_LIG.collapsed('time', iris.analysis.STD_DEV)
print mean_PI.data, mean_LIG.data, stdev_PI.data, stdev_LIG.data

months=SAM_PI.shape[0]
times= pd.date_range(start='1850-01-01', periods=months, freq='AS') #for annual values
#times= pd.date_range(start='1850-01-01', periods=months, freq='MS') #for monthly values
##edit times for select interval
#times=times[:50*12]
times=np.arange(0,months,1) #only number of show years from 1850 9for annual SAM)

zero_line=np.empty(SAM_PI.shape[0]); zero_line.fill(0)

fig = plt.figure()

ax = fig.add_subplot(211)

fig.autofmt_xdate()
plt.plot(times, zero_line, c='black' )
plt.plot(times, SAM_PI.data, c='black',linewidth=1, label='PI')
plt.fill_between(times, 0, np.ma.masked_where(SAM_PI.data <= 0, SAM_PI.data) , alpha=0.5, facecolor='orange')
plt.fill_between(times, 0, np.ma.masked_where(SAM_PI.data >= 0, SAM_PI.data) , alpha=0.5, facecolor='blue')
plt.xlabel('Time(years)', fontsize=14)
plt.ylabel('SAM Index', fontsize=14)
plt.title(r'SAM Index PI 1850-2050', fontsize=16)
#plt.ticklabel_format(useOffset=False)
#plt.ticklabel_format(axis='y',useOffset=False)  #only on y axis if we are using a date array created with pandas
#ax.xaxis.set_major_locator(mdates.YearLocator(20))
plt.legend()


ax = fig.add_subplot(212)

fig.autofmt_xdate()
plt.plot(times, zero_line, c='black' )
plt.plot(times, SAM_LIG.data, c='black',linewidth=1,  label='LIG')
plt.fill_between(times, 0, np.ma.masked_where(SAM_LIG.data <= 0, SAM_LIG.data) , alpha=0.5, facecolor='orange')
plt.fill_between(times, 0, np.ma.masked_where(SAM_LIG.data >= 0, SAM_LIG.data) , alpha=0.5, facecolor='blue')
plt.xlabel('Time(years)', fontsize=14)
plt.ylabel('SAM Index', fontsize=14)
plt.title(r'SAM Index LIG 1850-2050', fontsize=16)
#plt.ticklabel_format(useOffset=False)
#plt.ticklabel_format(axis='y',useOffset=False)  #only on y axis if we are using a date array created with pandas
#ax.xaxis.set_major_locator(mdates.YearLocator(20))
plt.legend()

#iplt.show()



###################### PLOT U at 850 hPa for SH JET POSITION AND STRENGTH ######################################

U850_PI=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/uar766_U_V_3D_18502050.nc', iris.Constraint('x_wind', pressure=lambda p: p == 300)) #850))      
U850_LIG=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/uba937_U_V_W_3D_18502050.nc', iris.Constraint('x_wind', pressure=lambda p: p == 300)) #850))
V850_LIG=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/uba937_U_V_W_3D_18502050.nc', iris.Constraint('y_wind', pressure=lambda p: p == 300)) #850))
#W_LIG=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/uba937_U_V_W_3D_18502050.nc', iris.Constraint('upward_air_velocity', pressure=lambda p: p == 850))

#computing 200yr means for PI and LIG
U850_LIG_mean=U850_LIG.collapsed('time', iris.analysis.MEAN)
V850_LIG_mean=V850_LIG.collapsed('time', iris.analysis.MEAN)
U850_PI_mean=U850_PI.collapsed('time', iris.analysis.MEAN)
#W_LIG_mean=W_LIG.collapsed('time', iris.analysis.MEAN)
#computing mean of LIG-PI diff and STdev of LIG for S/N
U850_mean_diff= (U850_LIG[:,:,:]-U850_PI[:,:,:]).collapsed('time', iris.analysis.MEAN)
U850_PI_stdev=U850_PI.collapsed('time', iris.analysis.STD_DEV)


fig=plt.figure(figsize=(20,5), tight_layout=True)
fig.suptitle('U at 850 hPa 200-year mean (m/s)' , x=0.5, y=0.05, fontsize=28)

ax = fig.add_subplot(141)
ax=plt.gca()
bmap=Basemap(projection='ortho',lat_0=-90,lon_0=0,resolution='l', ax=ax)

#get the lon,lat values from cube in order to plot: 
#lon=U850_LIG_mean.coord('longitude').points
lon=np.linspace(0, 360, 192) 
lat=U850_LIG_mean.coord('latitude').points
x,y=bmap(*np.meshgrid(lon,lat))

orig_cmap=matplotlib.cm.RdBu_r 
shifted_cmap=shiftedColorMap(orig_cmap, midpoint=0.40, name='shifted')
levels=np.arange(-10, 17, 1)
contours=bmap.contourf(x,y, U850_LIG_mean[:,:].data, levels, cmap='shifted') #cmap='RdBu_r') 

bmap.drawcoastlines()
bmap.drawparallels(np.arange(-80.,81.,20.))
bmap.drawmeridians(np.arange(-180.,181.,20.))
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label('U at 850 hPa', fontsize=16) 
plt.title(r'LIG', fontsize=16)

ax = fig.add_subplot(142)
ax=plt.gca()
bmap=Basemap(projection='ortho',lat_0=-90,lon_0=0,resolution='l', ax=ax)
#map.shadedrelief(scale=0.25)

#get the lon,lat values from cube in order to plot: 
lon=np.linspace(0, 360, 192) 
lat=U850_PI_mean.coord('latitude').points
x,y=bmap(*np.meshgrid(lon,lat))

orig_cmap=matplotlib.cm.RdBu_r 
shifted_cmap=shiftedColorMap(orig_cmap, midpoint=0.40, name='shifted')
contours=bmap.contourf(x,y, U850_PI_mean[:,:].data, levels, cmap='shifted') #cmap='RdBu_r') 

bmap.drawcoastlines()
bmap.drawparallels(np.arange(-80.,81.,20.))
bmap.drawmeridians(np.arange(-180.,181.,20.))
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label('U at 850 hPa', fontsize=16)
plt.title(r'PI', fontsize=16)


ax = fig.add_subplot(143)
ax=plt.gca()
bmap=Basemap(projection='ortho',lat_0=-90,lon_0=0,resolution='l', ax=ax)
#map.shadedrelief(scale=0.25)

orig_cmap=matplotlib.cm.PuOr_r
shifted_cmap=shiftedColorMap(orig_cmap, midpoint=0.58, name='shifted')
x,y=bmap(*np.meshgrid(lon,lat))
levels_2=np.arange(-2.5, 1.8, 0.2)
contours=bmap.contourf(x,y, (U850_mean_diff.data ) , levels_2, cmap='shifted') #cmap='PuOr_r') 
bmap.quiver(x, y, U850_LIG_mean[:,:].data,V850_LIG_mean[:,:].data, pivot='middle', units='xy') #, scale=1, scale_units='xy')

bmap.drawcoastlines()
bmap.drawparallels(np.arange(-80.,81.,20.))
bmap.drawmeridians(np.arange(-180.,181.,20.))
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label('U$_{LIG}$ - U$_{PI}$ at 850 hPa', fontsize=16) 
plt.title(r'LIG - PI', fontsize=16)

ax = fig.add_subplot(144)
ax=plt.gca()
bmap=Basemap(projection='ortho',lat_0=-90,lon_0=0,resolution='l', ax=ax)
#map.shadedrelief(scale=0.25)
x,y=bmap(*np.meshgrid(lon,lat))

#U850_PI_stdev.data=np.ma.masked_less_equal(U850_PI_stdev.data, 0.01)
S_to_N=U850_mean_diff/U850_PI_stdev
levels=np.arange(0, 1.1, 0.1)
contours=bmap.contourf(x,y, abs(S_to_N.data) , levels, extend='max', cmap='YlGnBu') 
#contours=bmap.contourf(x,y, U850_LIG_stdev.data , levels, cmap='GnBu_r') 

bmap.drawcoastlines()
bmap.drawparallels(np.arange(-80.,81.,20.))
bmap.drawmeridians(np.arange(-180.,181.,20.))
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label('S/N of U$_{LIG}$ - U$_{PI}$ at 850 hPa', fontsize=16) 
plt.title(r'S/N of LIG - PI', fontsize=16)

#plt.show()

###################### PLOT zonal mean cross section of U vs height/pressure  ######################################

U_PI=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/uar766_U_V_3D_18502050.nc', iris.Constraint('x_wind', latitude= lambda lat: -90 <= lat <= 0))       
U_LIG=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/uba937_U_V_W_3D_18502050.nc', iris.Constraint('x_wind', latitude= lambda lat: -90 <= lat <= 0))

U_LIG_zmn=U_LIG.collapsed('longitude', iris.analysis.MEAN)
U_PI_zmn=U_PI.collapsed('longitude', iris.analysis.MEAN)

U_LIG_mean=U_LIG_zmn.collapsed('time', iris.analysis.MEAN)
U_PI_mean=U_PI_zmn.collapsed('time', iris.analysis.MEAN)
#print U_LIG_mean

#computing mean of LIG-PI diff and STdev of LIG for S/N
U_mean_diff= (U_LIG_zmn[:,:,:]-U_PI_zmn[:,:,:]).collapsed('time', iris.analysis.MEAN)
U_PI_stdev=U_PI_zmn.collapsed('time', iris.analysis.STD_DEV)


fig=plt.figure(figsize=(20,5), tight_layout=True)
fig.suptitle(r'Zonal Mean Zonal Wind $\bar{U}$ 200-year mean (m/s)' , x=0.5, y=0.05, fontsize=28)


ax = fig.add_subplot(141)
ax=plt.gca()
contours=iplt.contourf(U_LIG_mean[:,:],cmap='RdBu_r') #,levels)
#contours_2=iplt.contour(U_LIG_mean[:,:], levels, colors='black', linewidths=1)
#plt.clabel(contours_2)
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label(r'$\bar{U}_{LIG}$', fontsize=16) 
plt.title(r'LIG', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
plt.ylabel('Height (hPa)', fontsize=14)
plt.xticks([-80, -60, -40, -20, 0])


ax = fig.add_subplot(142)
ax=plt.gca()
contours=iplt.contourf(U_PI_mean[:,:],cmap='RdBu_r') #,levels
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label(r'$\bar{U}_{PI}$', fontsize=16) 
plt.title(r'PI', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
#plt.ylabel('Height (hPa)', fontsize=14)
plt.xticks([-80, -60, -40, -20, 0])

ax = fig.add_subplot(143)
ax=plt.gca()
contours=iplt.contourf(U_mean_diff ,cmap='PuOr_r') #,levels
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label(r'$\bar{U}_{LIG} - \bar{U}_{PI}$', fontsize=16) 
plt.title(r'LIG - PI ', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
#plt.ylabel('Height (hPa)', fontsize=14)
plt.xticks([-80, -60, -40, -20, 0])


ax = fig.add_subplot(144)
ax=plt.gca()

S_to_N=(U_mean_diff/U_PI_stdev)
S_to_N.data=abs(S_to_N.data)
levels=np.arange(0, 1.1, 0.1)
contours=iplt.contourf(S_to_N , levels, extend='max', cmap='YlGnBu')
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
#cbar.set_label(r'$S/N of \bar{U}_{LIG} - \bar{U}_{PI}$', fontsize=16) 
plt.title(r'S/N of LIG - PI ', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
#plt.ylabel('Height (hPa)', fontsize=14)
plt.xticks([-80, -60, -40, -20, 0])



#plt.show()

###################### PLOT cross section of W  and U vs height/pressure  ######################################

W_LIG=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/uba937_U_V_W_3D_18502050.nc', iris.Constraint('upward_air_velocity', latitude= lambda lat: -90 <= lat <= 0))
W_PI=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/uas245_W_3D_18502050.nc', iris.Constraint('upward_air_velocity', latitude= lambda lat: -90 <= lat <= 0)) #NEED MO PI!!

V_PI=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/uar766_U_V_3D_18502050.nc', iris.Constraint('y_wind', latitude= lambda lat: -90 <= lat <= 0))       
V_LIG=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/uba937_U_V_W_3D_18502050.nc', iris.Constraint('y_wind', latitude= lambda lat: -90 <= lat <= 0))

W_LIG_mean=W_LIG.collapsed('time', iris.analysis.MEAN)
W_PI_mean=W_PI.collapsed('time', iris.analysis.MEAN)
cross=iris.Constraint(latitude=lambda lat: lat == -70, longitude= lambda lon: 0<= lon <= 359.125) #lat == -65,
W_LIG_mean=W_LIG_mean.extract(cross)
W_PI_mean=W_PI_mean.extract(cross)

U_LIG_mean=U_LIG.collapsed('time', iris.analysis.MEAN)
U_PI_mean=U_PI.collapsed('time', iris.analysis.MEAN)
cross=iris.Constraint(latitude=lambda lat: lat == -70, longitude= lambda lon: 0<= lon <= 359.125)
U_LIG_mean=U_LIG_mean.extract(cross)
U_PI_mean=U_PI_mean.extract(cross)
#print U_PI_mean

V_LIG_mean=V_LIG.collapsed('time', iris.analysis.MEAN)
V_PI_mean=V_PI.collapsed('time', iris.analysis.MEAN)
cross=iris.Constraint(latitude=lambda lat: lat == -70, longitude= lambda lon: 0<= lon <= 359.125)
V_LIG_mean=V_LIG_mean.extract(cross)
V_PI_mean=V_PI_mean.extract(cross)

#compute wind speed
#UV_LIG_mean=(U_LIG_mean**2 +V_LIG_mean**2)**0.5
#UV_PI_mean= (U_PI_mean**2 +V_PI_mean**2)**0.5

##preparing data to plot wind vectors
x=U_LIG_mean.coord('longitude').points 
y=U_LIG_mean.coord('pressure').points
##normalize wind vectors for uniform arrow size
U_norm_LIG= U_LIG_mean/(U_LIG_mean**2 + W_LIG_mean**2)**0.5
W_norm_LIG= W_LIG_mean/(U_LIG_mean**2 + W_LIG_mean**2)**0.5
V_norm_LIG= V_LIG_mean/(V_LIG_mean**2 + V_LIG_mean**2)**0.5
U_norm_PI= U_PI_mean/(U_PI_mean**2 + W_PI_mean**2)**0.5
W_norm_PI= W_PI_mean/(U_PI_mean**2 + W_PI_mean**2)**0.5
V_norm_PI= V_PI_mean/(V_PI_mean**2 + V_PI_mean**2)**0.5


fig=plt.figure(figsize=(20,5), tight_layout=True)
fig.suptitle(r'Zonal cross-section of W 200-year mean (m/s)' , x=0.5, y=0.05, fontsize=28)


ax = fig.add_subplot(231)
ax=plt.gca()
contours=iplt.contourf(W_LIG_mean[:,:],cmap='RdBu_r') #,levels
plt.quiver(x,y, U_norm_LIG[:,:].data, W_norm_LIG[:,:].data, pivot='middle', units='xy' ) #, scale=1, scale_units='xy')

ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
plt.title(r'W LIG', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
plt.ylabel('Height (hPa)', fontsize=14)
#plt.xticks([-80, -60, -40, -20, 0])

ax = fig.add_subplot(232)
ax=plt.gca()
contours=iplt.contourf(W_PI_mean[:,:],cmap='RdBu_r') #,levels)
plt.quiver(x,y, U_norm_PI[:,:].data, W_norm_PI[:,:].data, pivot='middle', units='xy') #, scale=1, scale_units='xy') 

ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
plt.title(r'W PI', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
plt.ylabel('Height (hPa)', fontsize=14)
#plt.xticks([-80, -60, -40, -20, 0])

ax = fig.add_subplot(233)
ax=plt.gca()
contours=iplt.contourf(W_LIG_mean[:,:]-W_PI_mean[:,:],cmap='seismic') #,levels)
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
plt.title(r' W LIG - W PI', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
plt.ylabel('Height (hPa)', fontsize=14)
#plt.xticks([-80, -60, -40, -20, 0])

ax = fig.add_subplot(234)
ax=plt.gca()
contours=iplt.contourf(U_LIG_mean[:,:],cmap='PuOr_r') #,levels)
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
plt.title(r'U LIG', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
plt.ylabel('Height (hPa)', fontsize=14)
#plt.xticks([-80, -60, -40, -20, 0])

ax = fig.add_subplot(235)
ax=plt.gca()
contours=iplt.contourf(U_PI_mean[:,:],cmap='PuOr_r') #,levels)
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
plt.title(r'U PI', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
plt.ylabel('Height (hPa)', fontsize=14)
#plt.xticks([-80, -60, -40, -20, 0])

ax = fig.add_subplot(236)
ax=plt.gca()
contours=iplt.contourf(U_LIG_mean[:,:]-U_PI_mean[:,:],cmap='seismic') #,levels)
ax.invert_yaxis()
cbar = plt.colorbar(contours, orientation='horizontal')
plt.title(r'U LIG - U PI', fontsize=16)
plt.xlabel('Latitude (Degrees)', fontsize=14)
plt.ylabel('Height (hPa)', fontsize=14)
#plt.xticks([-80, -60, -40, -20, 0])

plt.show()



