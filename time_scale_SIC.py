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

from shiftedColorMap import*
from scipy import stats



################## Define function to compute NH and SH sea ice area ##############################################################################################

def N_S(MO, AR):

 #Global average: 
 MO.coord('latitude').guess_bounds()
 MO.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(MO)
 seaice_area_AR=AR*cell_area
 seaice_area_MO=MO*cell_area


 seaice_global_MO= seaice_area_MO.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)

 seaice_global_AR= seaice_area_AR.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)


 #need to reset the coordinates bounds before to proceed
 MO.coord('latitude').bounds= None
 MO.coord('longitude').bounds= None

 #Extract SH
 SH=iris.Constraint(latitude=lambda lat: -89.5 <= lat <= 0, longitude= lambda lon: -179.5 <= lon <= 179.5)
 MO_SH=MO.extract(SH)
 AR_SH=AR.extract(SH)


 #spatial average:
 MO_SH.coord('latitude').guess_bounds()
 MO_SH.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(MO_SH)
 #sea ice area (aggregate) in input is aice (the sea ice area fraction), to calculate the sea ice area in m^2 we use the cell area: 
 seaice_area_AR=AR_SH*cell_area
 seaice_area_MO=MO_SH*cell_area


 seaice_SH_MO= seaice_area_MO.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)

 seaice_SH_AR= seaice_area_AR.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)

 #print seaice_SH_MO
 #need to reset the coordinates bounds before to proceed
 MO.coord('latitude').bounds= None
 MO.coord('longitude').bounds= None

 #Extract NH
 NH=iris.Constraint(latitude=lambda lat: 0 <= lat <= 89.5, longitude= lambda lon: -179.5 <= lon <= 179.5)
 MO_NH=MO.extract(NH)
 AR_NH=AR.extract(NH)

 #spatial average:
 MO_NH.coord('latitude').guess_bounds()
 MO_NH.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(MO_NH)
 #sea ice area (aggregate) in input is aice (the sea ice area fraction), to calculate the sea ice area in m^2 we use the cell area: 
 seaice_area_AR=AR_NH*cell_area
 seaice_area_MO=MO_NH*cell_area


 seaice_NH_MO= seaice_area_MO.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)

 seaice_NH_AR= seaice_area_AR.collapsed(['latitude', 'longitude'],
                                        iris.analysis.SUM)

 return seaice_global_MO, seaice_global_AR, seaice_NH_MO, seaice_NH_AR, seaice_SH_MO, seaice_SH_AR

################################################ define function to compute seasonal cycle #################################################################################################

def s_cycle(area_MO, area_AR):
 month=area_MO.shape[0]
 mean_MO_cl=iris.cube.CubeList()
 mean_AR_cl=iris.cube.CubeList()
 stdev_MO_cl=iris.cube.CubeList()
 stdev_AR_cl=iris.cube.CubeList()

 for i in range (0, 12): # in python's range (start , end) 'end' is not included (this is effectively: i=0,11)
   mean_MO_cl.append(area_MO[i:month:12].collapsed('time', iris.analysis.MEAN))
   mean_AR_cl.append(area_AR[i:month:12].collapsed('time', iris.analysis.MEAN))
   stdev_MO_cl.append(area_MO[i:month:12].collapsed('time', iris.analysis.STD_DEV))
   stdev_AR_cl.append(area_AR[i:month:12].collapsed('time', iris.analysis.STD_DEV))
     
 mean_MO= mean_MO_cl.merge_cube()
 mean_AR= mean_AR_cl.merge_cube()                   
 stdev_MO=stdev_MO_cl.merge_cube()
 stdev_AR=stdev_AR_cl.merge_cube()

 return mean_MO, mean_AR, stdev_MO, stdev_AR

###############################################################################################################################################################################

#load annual values
MO=iris.load_cube('/nerc/n02/n02/vittoria/seaice_annual_uar766_*.nc',  'ice area  (aggregate)' ) 
AR_cubes=iris.load('/nerc/n02/n02/vittoria/seaice_annual_uas245_*.nc',     'ice area  (aggregate)')
#concatenate cubes in one single dataset to plot
iris.util.unify_time_units(AR_cubes)
AR=AR_cubes.concatenate_cube()

#Now use N_S function on yearly data:
seaice_global_MO, seaice_global_AR, seaice_NH_MO, seaice_NH_AR, seaice_SH_MO, seaice_SH_AR = N_S(MO, AR)


##10 years differences 
mean_Gl_10yr=iris.cube.CubeList()
mean_NH_10yr=iris.cube.CubeList()
mean_SH_10yr=iris.cube.CubeList()
stdev_Gl_10yr=iris.cube.CubeList()
stdev_NH_10yr=iris.cube.CubeList()
stdev_SH_10yr=iris.cube.CubeList()

for i in range (10, 200+1, 10):
  st=i-10
  mean_Gl_10yr.append( (seaice_global_MO[st:i] - seaice_global_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_10yr.append( (seaice_global_MO[st:i] - seaice_global_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_10yr.append( (seaice_NH_MO[st:i] - seaice_NH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_10yr.append( (seaice_NH_MO[st:i] - seaice_NH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_10yr.append( (seaice_SH_MO[st:i] - seaice_SH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_10yr.append( (seaice_SH_MO[st:i] - seaice_SH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )

##30 years differences 
mean_Gl_30yr=iris.cube.CubeList()
mean_NH_30yr=iris.cube.CubeList()
mean_SH_30yr=iris.cube.CubeList()
stdev_Gl_30yr=iris.cube.CubeList()
stdev_NH_30yr=iris.cube.CubeList()
stdev_SH_30yr=iris.cube.CubeList()

for i in range (30, 200+1, 30):
  st=i-30
  mean_Gl_30yr.append( (seaice_global_MO[st:i] - seaice_global_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_30yr.append( (seaice_global_MO[st:i] - seaice_global_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_30yr.append( (seaice_NH_MO[st:i] - seaice_NH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_30yr.append( (seaice_NH_MO[st:i] - seaice_NH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_30yr.append( (seaice_SH_MO[st:i] - seaice_SH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_30yr.append( (seaice_SH_MO[st:i] - seaice_SH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )

##50 years differences 
mean_Gl_50yr=iris.cube.CubeList()
mean_NH_50yr=iris.cube.CubeList()
mean_SH_50yr=iris.cube.CubeList()
stdev_Gl_50yr=iris.cube.CubeList()
stdev_NH_50yr=iris.cube.CubeList()
stdev_SH_50yr=iris.cube.CubeList()

for i in range (50, 200+1, 50):
  st=i-50
  mean_Gl_50yr.append( (seaice_global_MO[st:i] - seaice_global_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_50yr.append( (seaice_global_MO[st:i] - seaice_global_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_50yr.append( (seaice_NH_MO[st:i] - seaice_NH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_50yr.append( (seaice_NH_MO[st:i] - seaice_NH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_50yr.append( (seaice_SH_MO[st:i] - seaice_SH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_50yr.append( (seaice_SH_MO[st:i] - seaice_SH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )


##100 years differences 
mean_Gl_100yr=iris.cube.CubeList()
mean_NH_100yr=iris.cube.CubeList()
mean_SH_100yr=iris.cube.CubeList()
stdev_Gl_100yr=iris.cube.CubeList()
stdev_NH_100yr=iris.cube.CubeList()
stdev_SH_100yr=iris.cube.CubeList()

for i in range (100, 200+1, 100):
  st=i-100
  mean_Gl_100yr.append( (seaice_global_MO[st:i] - seaice_global_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_Gl_100yr.append( (seaice_global_MO[st:i] - seaice_global_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_NH_100yr.append( (seaice_NH_MO[st:i] - seaice_NH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_NH_100yr.append( (seaice_NH_MO[st:i] - seaice_NH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  mean_SH_100yr.append( (seaice_SH_MO[st:i] - seaice_SH_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
  stdev_SH_100yr.append( (seaice_SH_MO[st:i] - seaice_SH_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )


##200 years differences 
mean_Gl_200yr=  (seaice_global_MO - seaice_global_AR).collapsed('time', iris.analysis.MEAN) 
stdev_Gl_200yr= (seaice_global_MO - seaice_global_AR).collapsed('time', iris.analysis.STD_DEV) 
mean_NH_200yr=  (seaice_NH_MO - seaice_NH_AR).collapsed('time', iris.analysis.MEAN) 
stdev_NH_200yr= (seaice_NH_MO - seaice_NH_AR).collapsed('time', iris.analysis.STD_DEV) 
mean_SH_200yr=  (seaice_SH_MO - seaice_SH_AR).collapsed('time', iris.analysis.MEAN) 
stdev_SH_200yr= (seaice_SH_MO - seaice_SH_AR).collapsed('time', iris.analysis.STD_DEV) 

mean_Global=iris.cube.CubeList([mean_Gl_10yr.merge_cube(), mean_Gl_30yr.merge_cube(), mean_Gl_50yr.merge_cube(), mean_Gl_100yr.merge_cube(), mean_Gl_200yr])
stdev_Global=iris.cube.CubeList([stdev_Gl_10yr.merge_cube(), stdev_Gl_30yr.merge_cube(), stdev_Gl_50yr.merge_cube(), stdev_Gl_100yr.merge_cube(), stdev_Gl_200yr])
mean_NH=iris.cube.CubeList([mean_NH_10yr.merge_cube(), mean_NH_30yr.merge_cube(), mean_NH_50yr.merge_cube(), mean_NH_100yr.merge_cube(), mean_NH_200yr])
stdev_NH=iris.cube.CubeList([stdev_NH_10yr.merge_cube(), stdev_NH_30yr.merge_cube(), stdev_NH_50yr.merge_cube(), stdev_NH_100yr.merge_cube(), stdev_NH_200yr])
mean_SH=iris.cube.CubeList([mean_SH_10yr.merge_cube(), mean_SH_30yr.merge_cube(), mean_SH_50yr.merge_cube(), mean_SH_100yr.merge_cube(), mean_SH_200yr])
stdev_SH=iris.cube.CubeList([stdev_SH_10yr.merge_cube(), stdev_SH_30yr.merge_cube(), stdev_SH_50yr.merge_cube(), stdev_SH_100yr.merge_cube(), stdev_SH_200yr])

#print seaice_global_MO.collapsed('time', iris.analysis.MEAN).data 
#print seaice_global_MO.collapsed('time', iris.analysis.STD_DEV).data 
#print seaice_global_AR.collapsed('time', iris.analysis.MEAN).data 
#print seaice_global_AR.collapsed('time', iris.analysis.STD_DEV).data 

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
years=MO.shape[0]
#times= pd.date_range(start='1850-01-01', periods=years, freq='AS') 
times=np.arange(0,years,1)

#plotting 
fig=plt.figure(figsize=(20,10))

##subplot1
ax = fig.add_subplot(221) 
ax=plt.gca() 
ax.set_title('NH Sea Ice Area')
#fig.autofmt_xdate()

plt.plot( times, (seaice_NH_MO.data)/10**6, c='grey',  linewidth=2, label='MONSOON' )
plt.plot( times, (seaice_NH_AR.data)/10**6, c='black', linestyle= '--', linewidth=2, label='ARCHER')
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

plt.plot( times, (seaice_SH_MO.data)/10**6, c='grey',  linewidth=2, label='MONSOON' )
plt.plot( times, (seaice_SH_AR.data)/10**6, c='black', linestyle= '--', linewidth=2, label='ARCHER')
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
plt.ylabel('SIA$_{MO}$ - SIA$_{AR}$ ($km^{2}$)', fontsize=14)
plt.legend((nh, sh), ( 'NH','SH'), loc=1, fontsize=13)
plt.xticks([0, 1, 2, 3, 4], [ 10,  30,  50,  100,  200]) 
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.add_artist(AnchoredText('d', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )


##subplot3

#load monthly data:
MO_cubes=iris.load('/nerc/n02/n02/vittoria/seaice_monthly_uar766_*.nc', 'ice area  (aggregate)' )
AR_cubes=iris.load('/nerc/n02/n02/vittoria/seaice_monthly_uas245_*.nc', 'ice area  (aggregate)')
#concatenate cubes in one single dataset to plot
iris.util.unify_time_units(AR_cubes)
iris.util.unify_time_units(MO_cubes)
AR=AR_cubes.concatenate_cube()
MO=MO_cubes.concatenate_cube()

#Use N_S function on monthly data:
seaice_global_MO, seaice_global_AR, seaice_NH_MO, seaice_NH_AR, seaice_SH_MO, seaice_SH_AR = N_S(MO, AR)
#now use s_cycle function to compute seasonal cycles over the NH and SH
mean_NH_MO, mean_NH_AR, stdev_NH_MO, stdev_NH_AR = s_cycle(seaice_NH_MO, seaice_NH_AR)
mean_SH_MO, mean_SH_AR, stdev_SH_MO, stdev_SH_AR = s_cycle(seaice_SH_MO, seaice_SH_AR)


ax = fig.add_subplot(223)
ax=plt.gca()
ax.set_title('Cycle')

x_ax=np.arange(0,12,1)
plt.errorbar(x_ax, (mean_NH_MO.data)/10**6, yerr=(stdev_NH_MO.data)/10**6, c='grey',  alpha=0.8, linewidth=2, label='MONSOON' )
plt.errorbar(x_ax, (mean_NH_AR.data)/10**6, yerr=(stdev_NH_AR.data)/10**6, c='black',  alpha=0.8,  linestyle= '--', linewidth=2, label='ARCHER') 
plt.legend(loc=4)

plt.errorbar(x_ax, (mean_SH_MO.data)/10**6, yerr=(stdev_SH_MO.data)/10**6, c='grey',  alpha=0.8, linewidth=2) #, label='MONSOON' )
plt.errorbar(x_ax, (mean_SH_AR.data)/10**6, yerr=(stdev_SH_AR.data)/10**6, c='black',  alpha=0.8,  linestyle= '--', linewidth=2) #, label='ARCHER') 

plt.title(r'NH and SH sea ice area seasonal cycle (200-year mean) ', fontsize=16)
plt.xlabel('Months', fontsize=14)
plt.ylabel('SIA ($km^{2}$)', fontsize=14)
plt.ticklabel_format(axis='y',useOffset=False)

plt.xticks([-1, 0,   1,  2,  3,  4,  5,  6,  7, 8, 9,  10, 11 ],['', 'J', 'F', 'M', 'A', 'M','J','J','A','S','O','N','D'])
ax.add_artist(AnchoredText('c', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )



###let's now plot the maximum values for positive and negative differences at each time scale in a log-log plot to evaluate the existence of a power law behaviour
fig2=plt.figure()

def power_law(mean, color, lb, mrk): 
 pos=([])
 neg=([])
 for n in range (0, 4): 
  mean_p=0.
  mean_n=0.
  for i in range(0,mean[n].shape[0]):
    #if i >= 1: #of positive values, let's take the maximum value for each time scale  
     if mean[n][i].data > 0 and (mean[n][i].data > mean_p) : 
       mean_p = mean[n][i].data
     #now of negative values, let's take the maximum value for each time scale  
     if mean[n][i].data < 0 and abs(mean[n][i].data) > abs(mean_n) : 
       mean_n = abs(mean[n][i].data)
  pos += [mean_p]
  neg += [mean_n]

 if mean[4].data > 0:
  pos += [mean[4].data]
  t_p=np.array([10, 30, 50, 100, 200])
  t_n=np.array([10, 30, 50, 100])
  #plot1=plt.plot(t_p, pos, c=color, marker= mrk, markersize=5, label=lb)
  plot1=plt.scatter(np.log10(t_p), np.log10(pos), c=color, marker= mrk, s=60, label=lb) 
  slope, intercept, r_value, p_value, std_err=stats.linregress(np.log10(t_p), np.log10(pos)) #let's fit a line to the log-log plot
  line=slope*np.log10(t_p) +intercept
  plt.plot(np.log10(t_p), line, c=color)
  print 'pos', pos 
  print 'slope', slope, r_value, std_err
 else:
  neg += [abs(mean[4].data)]
  t_n=np.array([10, 30, 50, 100, 200])
  t_p=np.array([10, 30, 50, 100])
  #plot2=plt.plot(t_n, neg, c=color, marker= mrk, markersize=5, label=lb )
  plot2=plt.scatter(np.log10(t_n), np.log10(neg), c=color, marker= mrk, s=60, label=lb ) 
  slope, intercept, r_value, p_value, std_err=stats.linregress(np.log10(t_n), np.log10(neg))
  line=slope*np.log10(t_n) +intercept
  plt.plot(np.log10(t_n), line, c=color)
  print 'neg', neg
  print 'slope', slope, r_value, std_err

 #print pos
 #print neg
##if we want to plot both positive and negative:
 #plot1=plt.plot(t_p, pos, c=color, marker= mrk, markersize=5, label=lb)
 #plot2=plt.plot(t_n, neg, marker= mrk, markersize=5, c=color )
 
 return

#power_law(mean_Global, color='green', lb = 'Global', mrk='o')
power_law(mean_NH, color='orange', lb = 'NH', mrk='D')
power_law(mean_SH, color='red', lb = 'SH', mrk='x')

plt.title(r'Max SIA differences as a function of the averaging timescale', fontsize=16)
plt.xlabel('Log of Timescale (years)', fontsize=14)
plt.ylabel('Log of SIA$_{MO}$ - SIA$_{AR}$ ($km^{2}$)', fontsize=14)
plt.legend(loc=1, fontsize=13)

#to be used with plt.plot to get a log log plot with no line fitting
#plt.xlim(2**-0.5, 2**2.5)
#plt.xscale('log', basex=2)
#plt.yscale('log', basey=2)


plt.show()
