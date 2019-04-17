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



################## Define function to compute NH and SH sea ice area ##############################################################################################

def N_S(PI, LIG):

 #Global average: 
 PI.coord('latitude').guess_bounds()
 PI.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(PI)
 seaice_area_LIG=LIG*cell_area
 seaice_area_PI=PI*cell_area


 seaice_global_PI= seaice_area_PI.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)

 seaice_global_LIG= seaice_area_LIG.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)


 #need to reset the coordinates bounds before to proceed
 PI.coord('latitude').bounds= None
 PI.coord('longitude').bounds= None

 #Extract SH
 SH=iris.Constraint(latitude=lambda lat: -89.5 <= lat <= 0, longitude= lambda lon: -179.5 <= lon <= 179.5)
 PI_SH=PI.extract(SH)
 LIG_SH=LIG.extract(SH)


 #spatial average:
 PI_SH.coord('latitude').guess_bounds()
 PI_SH.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(PI_SH)
 #sea ice area (aggregate) in input is aice (the sea ice area fraction), to calculate the sea ice area in m^2 we use the cell area: 
 seaice_area_LIG=LIG_SH*cell_area
 seaice_area_PI=PI_SH*cell_area


 seaice_SH_PI= seaice_area_PI.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)

 seaice_SH_LIG= seaice_area_LIG.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)

 #print seaice_SH_PI
 #need to reset the coordinates bounds before to proceed
 PI.coord('latitude').bounds= None
 PI.coord('longitude').bounds= None

 #Extract NH
 NH=iris.Constraint(latitude=lambda lat: 0 <= lat <= 89.5, longitude= lambda lon: -179.5 <= lon <= 179.5)
 PI_NH=PI.extract(NH)
 LIG_NH=LIG.extract(NH)

 #spatial average:
 PI_NH.coord('latitude').guess_bounds()
 PI_NH.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(PI_NH)
 #sea ice area (aggregate) in input is aice (the sea ice area fraction), to calculate the sea ice area in m^2 we use the cell area: 
 seaice_area_LIG=LIG_NH*cell_area
 seaice_area_PI=PI_NH*cell_area


 seaice_NH_PI= seaice_area_PI.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)

 seaice_NH_LIG= seaice_area_LIG.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)

 return seaice_global_PI, seaice_global_LIG, seaice_NH_PI, seaice_NH_LIG, seaice_SH_PI, seaice_SH_LIG

################################################ define function to compute seasonal cycle #################################################################################################

def s_cycle(area_PI, area_LIG):
 month=area_PI.shape[0]
 mean_PI_cl=iris.cube.CubeList()
 mean_LIG_cl=iris.cube.CubeList()
 stdev_PI_cl=iris.cube.CubeList()
 stdev_LIG_cl=iris.cube.CubeList()

 for i in range (0, 12): # in python's range (start , end) 'end' is not included (this is effectively: i=0,11)
   mean_PI_cl.append(area_PI[i:month:12].collapsed('time', iris.analysis.MEAN))
   mean_LIG_cl.append(area_LIG[i:month:12].collapsed('time', iris.analysis.MEAN))
   stdev_PI_cl.append(area_PI[i:month:12].collapsed('time', iris.analysis.STD_DEV))
   stdev_LIG_cl.append(area_LIG[i:month:12].collapsed('time', iris.analysis.STD_DEV))
     
 mean_PI= mean_PI_cl.merge_cube()
 mean_LIG= mean_LIG_cl.merge_cube()                   
 stdev_PI=stdev_PI_cl.merge_cube()
 stdev_LIG=stdev_LIG_cl.merge_cube()

 return mean_PI, mean_LIG, stdev_PI, stdev_LIG

###############################################################################################################################################################################

#load annual values
LIG_cubes=iris.load('/gws/nopw/j04/pmip4_vol1/users/vittoria/seaice_annual_uba937_*.nc',  'ice area  (aggregate)' ) 
PI=iris.load_cube('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/seaice_annual_uar766_18502050.nc',     'ice area  (aggregate)')
#concatenate cubes in one single dataset to plot
iris.util.unify_time_units(LIG_cubes)
LIG=LIG_cubes.concatenate_cube()

'''
##load Feb and Sep values
LIG=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/aice_feb_sep_uba937_18502050.nc',  'aice feb' ) # 'aice sep'
PI=iris.load_cube('/gws/nopw/j04/pmip4_vol1/users/vittoria/aice_feb_sep_uar766_18502050.nc',  'aice feb')
'''

#Now use N_S function on yearly data:
seaice_global_PI, seaice_global_LIG, seaice_NH_PI, seaice_NH_LIG, seaice_SH_PI, seaice_SH_LIG = N_S(PI, LIG)


##10 years differences 
mean_Gl_10yr=iris.cube.CubeList()
mean_NH_10yr=iris.cube.CubeList()
mean_SH_10yr=iris.cube.CubeList()
stdev_Gl_10yr=iris.cube.CubeList()
stdev_NH_10yr=iris.cube.CubeList()
stdev_SH_10yr=iris.cube.CubeList()

for i in range (10, 200+1, 10):
  st=i-10
  mean_Gl_10yr.append( (seaice_global_PI[st:i] - seaice_global_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_10yr.append( (seaice_global_PI[st:i] - seaice_global_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_10yr.append( (seaice_NH_PI[st:i] - seaice_NH_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_10yr.append( (seaice_NH_PI[st:i] - seaice_NH_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_10yr.append( (seaice_SH_PI[st:i] - seaice_SH_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_10yr.append( (seaice_SH_PI[st:i] - seaice_SH_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )

##30 years differences 
mean_Gl_30yr=iris.cube.CubeList()
mean_NH_30yr=iris.cube.CubeList()
mean_SH_30yr=iris.cube.CubeList()
stdev_Gl_30yr=iris.cube.CubeList()
stdev_NH_30yr=iris.cube.CubeList()
stdev_SH_30yr=iris.cube.CubeList()

for i in range (30, 200+1, 30):
  st=i-30
  mean_Gl_30yr.append( (seaice_global_PI[st:i] - seaice_global_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_30yr.append( (seaice_global_PI[st:i] - seaice_global_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_30yr.append( (seaice_NH_PI[st:i] - seaice_NH_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_30yr.append( (seaice_NH_PI[st:i] - seaice_NH_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_30yr.append( (seaice_SH_PI[st:i] - seaice_SH_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_30yr.append( (seaice_SH_PI[st:i] - seaice_SH_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )

##50 years differences 
mean_Gl_50yr=iris.cube.CubeList()
mean_NH_50yr=iris.cube.CubeList()
mean_SH_50yr=iris.cube.CubeList()
stdev_Gl_50yr=iris.cube.CubeList()
stdev_NH_50yr=iris.cube.CubeList()
stdev_SH_50yr=iris.cube.CubeList()

for i in range (50, 200+1, 50):
  st=i-50
  mean_Gl_50yr.append( (seaice_global_PI[st:i] - seaice_global_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_50yr.append( (seaice_global_PI[st:i] - seaice_global_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_50yr.append( (seaice_NH_PI[st:i] - seaice_NH_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_50yr.append( (seaice_NH_PI[st:i] - seaice_NH_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_50yr.append( (seaice_SH_PI[st:i] - seaice_SH_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_50yr.append( (seaice_SH_PI[st:i] - seaice_SH_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )


##100 years differences 
mean_Gl_100yr=iris.cube.CubeList()
mean_NH_100yr=iris.cube.CubeList()
mean_SH_100yr=iris.cube.CubeList()
stdev_Gl_100yr=iris.cube.CubeList()
stdev_NH_100yr=iris.cube.CubeList()
stdev_SH_100yr=iris.cube.CubeList()

for i in range (100, 200+1, 100):
  st=i-100
  mean_Gl_100yr.append( (seaice_global_PI[st:i] - seaice_global_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_100yr.append( (seaice_global_PI[st:i] - seaice_global_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_100yr.append( (seaice_NH_PI[st:i] - seaice_NH_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_100yr.append( (seaice_NH_PI[st:i] - seaice_NH_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_100yr.append( (seaice_SH_PI[st:i] - seaice_SH_LIG[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_100yr.append( (seaice_SH_PI[st:i] - seaice_SH_LIG[st:i]).collapsed('time', iris.analysis.STD_DEV) )


##200 years differences 
mean_Gl_200yr=  (seaice_global_PI - seaice_global_LIG).collapsed('time', iris.analysis.MEAN) 
stdev_Gl_200yr= (seaice_global_PI - seaice_global_LIG).collapsed('time', iris.analysis.STD_DEV) 
mean_NH_200yr=  (seaice_NH_PI - seaice_NH_LIG).collapsed('time', iris.analysis.MEAN) 
stdev_NH_200yr= (seaice_NH_PI - seaice_NH_LIG).collapsed('time', iris.analysis.STD_DEV) 
mean_SH_200yr=  (seaice_SH_PI - seaice_SH_LIG).collapsed('time', iris.analysis.MEAN) 
stdev_SH_200yr= (seaice_SH_PI - seaice_SH_LIG).collapsed('time', iris.analysis.STD_DEV) 

mean_Global=iris.cube.CubeList([mean_Gl_10yr.merge_cube(), mean_Gl_30yr.merge_cube(), mean_Gl_50yr.merge_cube(), mean_Gl_100yr.merge_cube(), mean_Gl_200yr])
stdev_Global=iris.cube.CubeList([stdev_Gl_10yr.merge_cube(), stdev_Gl_30yr.merge_cube(), stdev_Gl_50yr.merge_cube(), stdev_Gl_100yr.merge_cube(), stdev_Gl_200yr])
mean_NH=iris.cube.CubeList([mean_NH_10yr.merge_cube(), mean_NH_30yr.merge_cube(), mean_NH_50yr.merge_cube(), mean_NH_100yr.merge_cube(), mean_NH_200yr])
stdev_NH=iris.cube.CubeList([stdev_NH_10yr.merge_cube(), stdev_NH_30yr.merge_cube(), stdev_NH_50yr.merge_cube(), stdev_NH_100yr.merge_cube(), stdev_NH_200yr])
mean_SH=iris.cube.CubeList([mean_SH_10yr.merge_cube(), mean_SH_30yr.merge_cube(), mean_SH_50yr.merge_cube(), mean_SH_100yr.merge_cube(), mean_SH_200yr])
stdev_SH=iris.cube.CubeList([stdev_SH_10yr.merge_cube(), stdev_SH_30yr.merge_cube(), stdev_SH_50yr.merge_cube(), stdev_SH_100yr.merge_cube(), stdev_SH_200yr])

#print seaice_global_PI.collapsed('time', iris.analysis.MEAN).data 
#print seaice_global_PI.collapsed('time', iris.analysis.STD_DEV).data 
#print seaice_global_LIG.collapsed('time', iris.analysis.MEAN).data 
#print seaice_global_LIG.collapsed('time', iris.analysis.STD_DEV).data 

###### print total sea ice area for LIG and PI on 200 year mean and calculate % of reduction ########
PI_tot = seaice_SH_PI.collapsed('time', iris.analysis.MEAN)
LIG_tot = seaice_SH_LIG.collapsed('time', iris.analysis.MEAN)
reduction= (LIG_tot - PI_tot)/PI_tot
print PI_tot.data
print LIG_tot.data
print reduction.data
######################################################################################################


'''
#save in a csv file mean and stdev
fileout=open("stats_SIC.csv", "w") #'a'
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
years=PI.shape[0]
#times= pd.date_range(start='1850-01-01', periods=years, freq='AS') 
times=np.arange(0,years,1)

#plotting 
fig=plt.figure(figsize=(20,10))

##subplot1
ax = fig.add_subplot(221) 
ax=plt.gca() 
ax.set_title('NH Sea Ice Area')
#fig.autofmt_xdate()

plt.plot( times, (seaice_NH_PI.data)/10**6, c='aqua',  linewidth=3, label='PI' )
plt.plot( times, (seaice_NH_LIG.data)/10**6, c='orange', linestyle= '--', linewidth=3, label='LIG')
plt.title(r'Northern Hemisphere sea ice area', fontsize=16)
plt.xlabel('Time (years)', fontsize=14)
plt.ylabel('SIA ($km^{2}$)', fontsize=14)
plt.ticklabel_format(axis='y',useOffset=False)
#ax.xaxis.set_major_locator(mdates.YearLocator(20))
ax.add_artist(AnchoredText('a', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )
plt.legend(loc=2)


##subplot2
ax = fig.add_subplot(222)
ax=plt.gca()
ax.set_title('SH Sea Ice Area')
#fig.autofmt_xdate()

plt.plot( times, (seaice_SH_PI.data)/10**6, c='aqua',  linewidth=3, label='PI' )
plt.plot( times, (seaice_SH_LIG.data)/10**6, c='orange', linestyle= '--', linewidth=3, label='LIG')
plt.title(r'Southern Hemisphere sea ice area', fontsize=16)
plt.xlabel('Time (years)', fontsize=14)
plt.ylabel('SIA ($km^{2}$)', fontsize=14)
plt.ticklabel_format(axis='y',useOffset=False)
#ax.xaxis.set_major_locator(mdates.YearLocator(20))
ax.add_artist(AnchoredText('b', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )
plt.legend(loc=2)


##subplot4
ax = fig.add_subplot(224)
ax=plt.gca()
ax.set_title('Delta')
ax.axhline(y=0 , color='black', linestyle='--') #plot zero line

#t_scales=np.array([ [10]*20, [30]*6, [50]*4, [100]*2, [200] ])

for n in range (0, 4): 
 for i in range(0,mean_NH[n].shape[0]):
  nh=plt.scatter(n-0.08, (mean_NH[n][i].data)/10**6, marker='D',  s=60, facecolors='none', edgecolors='orange' )
  #plt.errorbar(n,  mean_NH[n][i].data, yerr= stdev_NH[n][i].data, color='orange', alpha=0.8, linewidth=0.5 )
plt.scatter(4-0.08, (mean_NH[4].data)/10**6, marker='D',  s=60, facecolors='none', edgecolors='orange'  )
#plt.errorbar(4,  mean_NH[4].data, yerr= stdev_NH[4].data, color='orange', alpha=0.8, linewidth=0.5 )


for n in range (0, 4): 
 for i in range(0,mean_SH[n].shape[0]):
  sh=plt.scatter(n+0.08, (mean_SH[n][i].data)/10**6, c='red' , marker='x',  s=60)
plt.scatter(4+0.08, (mean_SH[4].data)/10**6, c='red', marker='x',  s=60 )


plt.title(r'SIA differences as a function of the averaging timescale', fontsize=16)
plt.xlabel('Timescale (years)', fontsize=14)
plt.ylabel('SIA$_{PI}$ - SIA$_{LIG}$ ($km^{2}$)', fontsize=14)
plt.legend((nh, sh), ( 'NH','SH'), loc=1, fontsize=13)
plt.xticks([0, 1, 2, 3, 4], [ 10,  30,  50,  100,  200]) 
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.add_artist(AnchoredText('d', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )


##subplot3

#load monthly data:
PI_cubes=iris.load('/group_workspaces/jasmin4/bas_palaeoclim/vittoria/seaice_monthly_uar766_*.nc', 'ice area  (aggregate)' )
LIG_cubes=iris.load('/gws/nopw/j04/pmip4_vol1/users/vittoria/seaice_monthly_uba937_*.nc',  'ice area  (aggregate)')
#concatenate cubes in one single dataset to plot
iris.util.unify_time_units(LIG_cubes)
iris.util.unify_time_units(PI_cubes)
LIG=LIG_cubes.concatenate_cube()
PI=PI_cubes.concatenate_cube()

#Use N_S function on monthly data:
seaice_global_PI, seaice_global_LIG, seaice_NH_PI, seaice_NH_LIG, seaice_SH_PI, seaice_SH_LIG = N_S(PI, LIG)
#now use s_cycle function to compute seasonal cycles over the NH and SH
mean_NH_PI, mean_NH_LIG, stdev_NH_PI, stdev_NH_LIG = s_cycle(seaice_NH_PI, seaice_NH_LIG)
mean_SH_PI, mean_SH_LIG, stdev_SH_PI, stdev_SH_LIG = s_cycle(seaice_SH_PI, seaice_SH_LIG)


ax = fig.add_subplot(223)
ax=plt.gca()
ax.set_title('Cycle')

x_ax=np.arange(0,12,1)

plt.errorbar(x_ax, (mean_NH_PI.data)/10**6, yerr=(stdev_NH_PI.data)/10**6, c='acqua',  alpha=0.8, linewidth=2) #, label='PI' )
plt.errorbar(x_ax, (mean_NH_LIG.data)/10**6, yerr=(stdev_NH_LIG.data)/10**6, c='orange',  alpha=0.8,  linestyle= '--', linewidth=2) #, label='LIG') 

plt.errorbar(x_ax, (mean_SH_PI.data)/10**6, yerr=(stdev_SH_PI.data)/10**6, c='aqua',  alpha=0.8, linewidth=3,  label='PI') 
plt.errorbar(x_ax, (mean_SH_LIG.data)/10**6, yerr=(stdev_SH_LIG.data)/10**6, c='orange',  alpha=0.8,  linestyle= '--', linewidth=3, label='LIG') 
plt.legend(loc=4)

plt.title(r'NH and SH sea ice area seasonal cycle (200-year mean) ', fontsize=16)
#plt.title(r'SH sea ice area seasonal cycle (200-year mean) ', fontsize=16)
plt.xlabel('Months', fontsize=14)
plt.ylabel('SIA ($km^{2}$)', fontsize=14)
plt.ticklabel_format(axis='y',useOffset=False)

plt.xticks([-1, 0,   1,  2,  3,  4,  5,  6,  7, 8, 9,  10, 11 ],['', 'J', 'F', 'M', 'A', 'M','J','J','A','S','O','N','D'])
ax.add_artist(AnchoredText('c', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )


plt.show()
