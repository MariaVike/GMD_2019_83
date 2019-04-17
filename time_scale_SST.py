import iris 
import iris.plot as iplt
import iris.quickplot as qplt

import iris.analysis.cartography 
import cartopy.crs as ccrs
import cartopy.feature as cfe
from mpl_toolkits.basemap import Basemap, maskoceans

import numpy as np
#np.set_printoptions(threshold=np.inf)
import matplotlib
import matplotlib.pyplot as plt
from iris.experimental.equalise_cubes import *
from iris.experimental.animate import*
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
from matplotlib.offsetbox import AnchoredText

import warnings
import glob
import os 

import pandas as pd
import scipy

#from shiftedColorMap import*
from scipy import stats



########################## Define function to compute Global, NH and SH SST ################################################################################################

def G_N_S(MO, AR):
 #Global average: 

 MO.coord('latitude').guess_bounds()
 MO.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(MO)
 #print cell_area 
 #fileout=open("area_weights.txt", "w")
 #fileout.write("2D lat/lon area weights array ")
 #fileout.write(str(cell_area))


 sst_global_MO= MO.collapsed(['latitude', 'longitude'],
                                        iris.analysis.MEAN,
                                        weights=cell_area)

 sst_global_AR= AR.collapsed(['latitude', 'longitude'],
                                        iris.analysis.MEAN,
                                         weights=cell_area)


 #need to reset the coordinates bounds before to proceed
 MO.coord('latitude').bounds= None
 MO.coord('longitude').bounds= None

 #Extract SH
 SH=iris.Constraint(latitude=lambda lat: -89.5 <= lat <= 0, longitude= lambda lon: -179.5 <= lon <= 179.5)
 MO_SH=MO.extract(SH)
 AR_SH=AR.extract(SH)

 #qplt.contourf(MO_annual[0,:,:])
 #plt.show()

 #spatial average:
 MO_SH.coord('latitude').guess_bounds()
 MO_SH.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(MO_SH)


 sst_SH_MO= MO_SH.collapsed(['latitude', 'longitude'],
                                        iris.analysis.MEAN,
                                         weights=cell_area)

 sst_SH_AR= AR_SH.collapsed(['latitude', 'longitude'],
                                           iris.analysis.MEAN,
                                           weights=cell_area)

 #print sst_SH_MO

 #need to reset the coordinates bounds before to proceed
 MO.coord('latitude').bounds= None
 MO.coord('longitude').bounds= None

 #Extract NH
 NH=iris.Constraint(latitude=lambda lat: 0 <= lat <= 89.5, longitude= lambda lon: -179.5 <= lon <= 179.5)
 MO_NH=MO.extract(NH)
 AR_NH=AR.extract(NH)

 #qplt.contourf(MO_annual[0,:,:])
 #plt.show()

 #spatial average:
 MO_NH.coord('latitude').guess_bounds()
 MO_NH.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(MO_NH)


 sst_NH_MO= MO_NH.collapsed(['latitude', 'longitude'],
                                        iris.analysis.MEAN,
                                         weights=cell_area)

 sst_NH_AR= AR_NH.collapsed(['latitude', 'longitude'],
                                           iris.analysis.MEAN,
                                           weights=cell_area)

 #print sst_NH_MO

 return sst_global_MO, sst_global_AR, sst_NH_MO, sst_NH_AR, sst_SH_MO, sst_SH_AR

##########################################################################################################################################################################
#'''
#load annual SST data:
AR_cubes=iris.load('/gws/nopw/j04/pmip4_vol1/users/vittoria/seaice_annual_uba937_*.nc',  'sea surface temperature' ) #LIG
MO=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/seaice_annual_uar766_18502050.nc',     'sea surface temperature')#PI
#concatenate cubes in one single dataset to plot
iris.util.unify_time_units(AR_cubes)
AR=AR_cubes.concatenate_cube()
#'''
'''
#load Air temperature data 
MO=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/uar766_AirTemp_18502050.nc',  'air_temperature' ) #PI
AR=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/uba937_SAT15m_18502050.nc',     'air_temperature') #LIG
'''

#Now use G_N_S function on yearly data:
sst_global_MO, sst_global_AR, sst_NH_MO, sst_NH_AR, sst_SH_MO, sst_SH_AR = G_N_S(MO, AR)


##10 years differences 
mean_Gl_10yr=iris.cube.CubeList()
mean_NH_10yr=iris.cube.CubeList()
mean_SH_10yr=iris.cube.CubeList()
stdev_Gl_10yr=iris.cube.CubeList()
stdev_NH_10yr=iris.cube.CubeList()
stdev_SH_10yr=iris.cube.CubeList()

for i in range (10, 200+1, 10):
  st=i-10
  mean_Gl_10yr.append( (sst_global_MO[st:i] - sst_global_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_10yr.append( (sst_global_MO[st:i] - sst_global_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_10yr.append( (sst_NH_MO[st:i] - sst_NH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_10yr.append( (sst_NH_MO[st:i] - sst_NH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_10yr.append( (sst_SH_MO[st:i] - sst_SH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_10yr.append( (sst_SH_MO[st:i] - sst_SH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
#print mean_Gl_10yr
#print mean_Gl_10yr[0]
#print mean_Gl_10yr[-1]

##30 years differences 
mean_Gl_30yr=iris.cube.CubeList()
mean_NH_30yr=iris.cube.CubeList()
mean_SH_30yr=iris.cube.CubeList()
stdev_Gl_30yr=iris.cube.CubeList()
stdev_NH_30yr=iris.cube.CubeList()
stdev_SH_30yr=iris.cube.CubeList()

for i in range (30, 200+1, 30):
  st=i-30
  mean_Gl_30yr.append( (sst_global_MO[st:i] - sst_global_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_30yr.append( (sst_global_MO[st:i] - sst_global_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_30yr.append( (sst_NH_MO[st:i] - sst_NH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_30yr.append( (sst_NH_MO[st:i] - sst_NH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_30yr.append( (sst_SH_MO[st:i] - sst_SH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_30yr.append( (sst_SH_MO[st:i] - sst_SH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )

##50 years differences 
mean_Gl_50yr=iris.cube.CubeList()
mean_NH_50yr=iris.cube.CubeList()
mean_SH_50yr=iris.cube.CubeList()
stdev_Gl_50yr=iris.cube.CubeList()
stdev_NH_50yr=iris.cube.CubeList()
stdev_SH_50yr=iris.cube.CubeList()

for i in range (50, 200+1, 50):
  st=i-50
  mean_Gl_50yr.append( (sst_global_MO[st:i] - sst_global_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_50yr.append( (sst_global_MO[st:i] - sst_global_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_50yr.append( (sst_NH_MO[st:i] - sst_NH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_50yr.append( (sst_NH_MO[st:i] - sst_NH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_50yr.append( (sst_SH_MO[st:i] - sst_SH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_50yr.append( (sst_SH_MO[st:i] - sst_SH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )


##100 years differences 
mean_Gl_100yr=iris.cube.CubeList()
mean_NH_100yr=iris.cube.CubeList()
mean_SH_100yr=iris.cube.CubeList()
stdev_Gl_100yr=iris.cube.CubeList()
stdev_NH_100yr=iris.cube.CubeList()
stdev_SH_100yr=iris.cube.CubeList()

for i in range (100, 200+1, 100):
  st=i-100
  mean_Gl_100yr.append( (sst_global_MO[st:i] - sst_global_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_100yr.append( (sst_global_MO[st:i] - sst_global_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_100yr.append( (sst_NH_MO[st:i] - sst_NH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_100yr.append( (sst_NH_MO[st:i] - sst_NH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_100yr.append( (sst_SH_MO[st:i] - sst_SH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_100yr.append( (sst_SH_MO[st:i] - sst_SH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )


##200 years differences 
mean_Gl_200yr= (sst_global_MO - sst_global_AR).collapsed('time', iris.analysis.MEAN) 
stdev_Gl_200yr= (sst_global_MO - sst_global_AR).collapsed('time', iris.analysis.STD_DEV) 
mean_NH_200yr=  (sst_NH_MO - sst_NH_AR).collapsed('time', iris.analysis.MEAN) 
stdev_NH_200yr= (sst_NH_MO - sst_NH_AR).collapsed('time', iris.analysis.STD_DEV) 
mean_SH_200yr=  (sst_SH_MO - sst_SH_AR).collapsed('time', iris.analysis.MEAN) 
stdev_SH_200yr= (sst_SH_MO - sst_SH_AR).collapsed('time', iris.analysis.STD_DEV) 

mean_Global=iris.cube.CubeList([mean_Gl_10yr.merge_cube(), mean_Gl_30yr.merge_cube(), mean_Gl_50yr.merge_cube(), mean_Gl_100yr.merge_cube(), mean_Gl_200yr])
stdev_Global=iris.cube.CubeList([stdev_Gl_10yr.merge_cube(), stdev_Gl_30yr.merge_cube(), stdev_Gl_50yr.merge_cube(), stdev_Gl_100yr.merge_cube(), stdev_Gl_200yr])
mean_NH=iris.cube.CubeList([mean_NH_10yr.merge_cube(), mean_NH_30yr.merge_cube(), mean_NH_50yr.merge_cube(), mean_NH_100yr.merge_cube(), mean_NH_200yr])
stdev_NH=iris.cube.CubeList([stdev_NH_10yr.merge_cube(), stdev_NH_30yr.merge_cube(), stdev_NH_50yr.merge_cube(), stdev_NH_100yr.merge_cube(), stdev_NH_200yr])
mean_SH=iris.cube.CubeList([mean_SH_10yr.merge_cube(), mean_SH_30yr.merge_cube(), mean_SH_50yr.merge_cube(), mean_SH_100yr.merge_cube(), mean_SH_200yr])
stdev_SH=iris.cube.CubeList([stdev_SH_10yr.merge_cube(), stdev_SH_30yr.merge_cube(), stdev_SH_50yr.merge_cube(), stdev_SH_100yr.merge_cube(), stdev_SH_200yr])

#print sst_global_MO.collapsed('time', iris.analysis.MEAN).data 
#print sst_global_MO.collapsed('time', iris.analysis.STD_DEV).data 
#print sst_global_AR.collapsed('time', iris.analysis.MEAN).data 
#print sst_global_AR.collapsed('time', iris.analysis.STD_DEV).data 

'''
#save in a csv file mean and stdev
fileout=open("stats_SST.csv", "w") #'a'
#fileout=open("stats_Temp.csv", "w") #'a'
fileout.write("mean_10yr_Global" + "stdev_10yr_Global")
fileout.write(str(mean_Global[0].data) + str(stdev_Global[0].data))
fileout.write("mean_10yr_NH" + "stdev_10yr_NH")
fileout.write(str(mean_NH[0].data) + str(stdev_NH[0].data))
fileout.write("mean_10yr_SH" + "stdev_10yr_SH")
fileout.write(str(mean_SH[0].data) + str(stdev_SH[0].data))
fileout.write("mean_30yr_Global" + "stdev_30yr_Global")
fileout.write(str(mean_Global[1].data) + str(stdev_Global[1].data))
fileout.write("mean_30yr_NH" + "stdev_30yr_NH")
fileout.write(str(mean_NH[1].data) + str(stdev_NH[1].data))
fileout.write("mean_30yr_SH" + "stdev_30yr_SH")
fileout.write(str(mean_SH[1].data) + str(stdev_SH[1].data))
fileout.write("mean_50yr_Global" + "stdev_50yr_Global")
fileout.write(str(mean_Global[2].data) + str(stdev_Global[2].data))
fileout.write("mean_50yr_NH" + "stdev_50yr_NH")
fileout.write(str(mean_NH[2].data) + str(stdev_NH[2].data))
fileout.write("mean_50yr_SH" + "stdev_50yr_SH")
fileout.write(str(mean_SH[2].data) + str(stdev_SH[2].data))
fileout.write("mean_100yr_Global" + "stdev_100yr_Global")
fileout.write(str(mean_Global[3].data) + str(stdev_Global[3].data))
fileout.write("mean_100yr_NH" + "stdev_100yr_NH")
fileout.write(str(mean_NH[3].data) + str(stdev_NH[3].data))
fileout.write("mean_100yr_SH" + "stdev_100yr_SH")
fileout.write(str(mean_SH[3].data) + str(stdev_SH[3].data))
fileout.write("mean_200yr_Global" + "stdev_200yr_Global")
fileout.write(str(mean_Global[4].data) + str(stdev_Global[4].data))
fileout.write("mean_200yr_NH" + "stdev_200yr_NH")
fileout.write(str(mean_NH[4].data) + str(stdev_NH[4].data))
fileout.write("mean_200yr_SH" + "stdev_200yr_SH")
fileout.write(str(mean_SH[4].data) + str(stdev_SH[4].data))
'''

##create time coordinate for plots: 
years=MO.shape[0]
#times= pd.date_range(start='1850-01-01', periods=years, freq='AS') 
times=np.arange(0,years,1)

#plotting 
fig=plt.figure(figsize=(20,10))

#subplot1
ax = fig.add_subplot(221) 
ax=plt.gca()
ax.set_title('Global SST')
#fig.autofmt_xdate()

plt.plot( times, sst_global_MO.data, c='darkgreen',  linewidth=3, label='PI' )
plt.plot( times, sst_global_AR.data, c='orange', linestyle= '--', linewidth=3, label='LIG')
plt.title(r'Global sea surface temperature', fontsize=16)
#plt.title(r'Global air temperature', fontsize=16)
plt.xlabel('Time (years)', fontsize=14)
plt.ylabel('SST (DegC)', fontsize=14)
#plt.ylabel('SAT (DegC)', fontsize=14)
plt.ticklabel_format(axis='y',useOffset=False)
#ax.xaxis.set_major_locator(mdates.YearLocator(20))
ax.add_artist(AnchoredText('a', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )
plt.legend(loc=2)



##subplot2
ax = fig.add_subplot(222) 
ax=plt.gca() 
ax.set_title('NH SST')
#fig.autofmt_xdate()

plt.plot( times, sst_NH_MO.data, c='darkgreen',  linewidth=3, label='PI' )
plt.plot( times, sst_NH_AR.data, c='orange', linestyle= '--', linewidth=3, label='LIG')
plt.title(r'Northern Hemisphere sea surface temperature', fontsize=16)
#plt.title(r'Northern Hemisphere air temperature', fontsize=16)
plt.xlabel('Time (years)', fontsize=14)
plt.ylabel('SST (DegC)', fontsize=14)
#plt.ylabel('SAT (DegC)', fontsize=14)
plt.ticklabel_format(axis='y',useOffset=False)
#ax.xaxis.set_major_locator(mdates.YearLocator(20))
ax.add_artist(AnchoredText('b', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )
plt.legend(loc=2)


##subplot3
ax = fig.add_subplot(223)
ax=plt.gca()
ax.set_title('SH SST')
#fig.autofmt_xdate()

plt.plot( times, sst_SH_MO.data, c='darkgreen',  linewidth=3, label='PI' )
plt.plot( times, sst_SH_AR.data, c='orange', linestyle= '--', linewidth=3, label='LIG')
plt.title(r'Southern Hemisphere sea surface temperature', fontsize=16)
#plt.title(r'Southern Hemisphere air temperature', fontsize=16)
plt.xlabel('Time (years)', fontsize=14)
plt.ylabel('SST (DegC)', fontsize=14)
#plt.ylabel('SAT (DegC)', fontsize=14)
plt.ticklabel_format(axis='y',useOffset=False)
#ax.xaxis.set_major_locator(mdates.YearLocator(20))
ax.add_artist(AnchoredText('c', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )
plt.legend(loc=2)


##subplot4
ax = fig.add_subplot(224)
ax=plt.gca()
ax.set_title('Delta')
ax.axhline(y=0 , color='black', linestyle='--') #plot zero line

x_ax=np.array([10, 30, 50, 100, 200])

for n in range (0, 4): 
 for i in range(0,mean_NH[n].shape[0]):
  #nh=plt.scatter(x_ax[n]-5, mean_NH[n][i].data, marker='D',  s=60, facecolors='none', edgecolors='orange' ) #show real axis 
  nh=plt.scatter(n-0.15, mean_NH[n][i].data, marker='D',  s=60, facecolors='none', edgecolors='orange' ) # x_axis equally spaced
  #plt.errorbar(n,  mean_NH[n][i].data, yerr= stdev_NH[n][i].data, color='orange', alpha=0.8, linewidth=0.5 )

#the 200 yr averaging period is not a list, thus simply plot it by: 
#plt.scatter(x_ax[4]-5, mean_NH[4].data, marker='D',  s=60, facecolors='none', edgecolors='orange'  ) #show real axis 
plt.scatter(4-0.15, mean_NH[4].data, marker='D',  s=60, facecolors='none', edgecolors='orange'  )
#plt.errorbar(4,  mean_NH[4].data, yerr= stdev_NH[4].data, color='orange', alpha=0.8, linewidth=0.5 )

for n in range (0, 4): 
 for i in range(0,mean_Global[n].shape[0]):
  #gl=plt.scatter(x_ax[n], mean_Global[n][i].data, marker='o', s=60, facecolors='none', edgecolors='green')
   gl=plt.scatter([n], mean_Global[n][i].data, marker='o', s=60, facecolors='none', edgecolors='green')
#plt.scatter(x_ax[4], mean_Global[4].data, marker='o', s=60, facecolors='none', edgecolors='green')
plt.scatter(4, mean_Global[4].data, marker='o', s=60, facecolors='none', edgecolors='green')

for n in range (0, 4): 
 for i in range(0,mean_SH[n].shape[0]):
  #sh=plt.scatter(x_ax[n]+5, mean_SH[n][i].data, c='red' , marker='x',  s=60)
  sh=plt.scatter(n+0.15, mean_SH[n][i].data, c='red' , marker='x',  s=60)
#plt.scatter(x_ax[4]+5, mean_SH[4].data, c='red', marker='x',  s=60 )
plt.scatter(4+0.15, mean_SH[4].data, c='red', marker='x',  s=60 )

plt.title(r'SST differences as a function of the averaging timescale', fontsize=16)
#plt.title(r'SAT differences as a function of the averaging timescale', fontsize=16)
plt.xlabel('Timescale (years)', fontsize=14)
plt.ylabel('SST$_{MO}$ - SST$_{AR}$ (DegC)', fontsize=14)
#plt.ylabel('SAT$_{PI}$ - SAT$_{LIG}$ (DegC)', fontsize=14)
plt.legend((gl, nh, sh), ('Global', 'NH', 'SH'), loc=1, fontsize=13)
ax.add_artist(AnchoredText('d', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )
plt.xticks([0, 1, 2, 3, 4], [ 10,  30,  50,  100,  200]) 

plt.show()
