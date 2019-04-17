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



########################## Define function to compute Tropical and NH/SH Extratropical TOA LW and SW flux ################################################################################

def G_N_S(MO, AR):
 #Extract Tropics :

 Tr=iris.Constraint(latitude=lambda lat: -23.5 <= lat <= 23.5, longitude= lambda lon: -179.5 <= lon <= 179.5) #between Tropic of Cancer and Tropic of Capricorn
 MO_T=MO.extract(Tr)
 AR_T=AR.extract(Tr) 

 MO_T.coord('latitude').guess_bounds()
 MO_T.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(MO_T)


 tropics_MO= MO_T.collapsed(['latitude', 'longitude'],
                                        iris.analysis.MEAN,
                                        weights=cell_area)

 tropics_AR= AR_T.collapsed(['latitude', 'longitude'],
                                        iris.analysis.MEAN,
                                         weights=cell_area)


 #need to reset the coordinates bounds before to proceed
 MO.coord('latitude').bounds= None
 MO.coord('longitude').bounds= None


 #Extract SH Extratropics 
 SH_Ex=iris.Constraint(latitude=lambda lat: -23.5 < lat <= 66.5, longitude= lambda lon: -179.5 <= lon <= 179.5) #between Tropic of Capricorn and Antarctic Circle
 MO_SH_Ex=MO.extract(SH_Ex)
 AR_SH_Ex=AR.extract(SH_Ex)

 #spatial average:
 MO_SH_Ex.coord('latitude').guess_bounds()
 MO_SH_Ex.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(MO_SH_Ex)


 SH_Ex_MO= MO_SH_Ex.collapsed(['latitude', 'longitude'],
                                        iris.analysis.MEAN,
                                         weights=cell_area)

 SH_Ex_AR= AR_SH_Ex.collapsed(['latitude', 'longitude'],
                                           iris.analysis.MEAN,
                                           weights=cell_area)

 #print sst_SH_MO

 #need to reset the coordinates bounds before to proceed
 MO.coord('latitude').bounds= None
 MO.coord('longitude').bounds= None

 #Extract NH Extratropics 
 NH_Ex=iris.Constraint(latitude=lambda lat: 23.5 < lat <= 66.5, longitude= lambda lon: -179.5 <= lon <= 179.5) #between Tropic of Cancer and Arctic Circle
 MO_NH_Ex=MO.extract(NH_Ex)
 AR_NH_Ex=AR.extract(NH_Ex)

 #qplt.contourf(MO_annual[0,:,:])
 #plt.show()

 #spatial average:
 MO_NH_Ex.coord('latitude').guess_bounds()
 MO_NH_Ex.coord('longitude').guess_bounds()
 cell_area = iris.analysis.cartography.area_weights(MO_NH_Ex)


 NH_Ex_MO= MO_NH_Ex.collapsed(['latitude', 'longitude'],
                                        iris.analysis.MEAN,
                                         weights=cell_area)

 NH_Ex_AR= AR_NH_Ex.collapsed(['latitude', 'longitude'],
                                           iris.analysis.MEAN,
                                           weights=cell_area)

 #print sst_NH_MO

 return tropics_MO, tropics_AR, NH_Ex_MO, NH_Ex_AR, SH_Ex_MO, SH_Ex_AR

##########################################################################################################################################################################

########################## Define function to compute meand and stdev of MO - AR differences ################################################################################

def MO_minus_AR(tropics_MO, tropics_AR, NH_Ex_MO, NH_Ex_AR, SH_Ex_MO, SH_Ex_AR, fname="filename"):
 
 ##10 years differences 
 mean_Tr_10yr=iris.cube.CubeList()
 mean_NH_Ex_10yr=iris.cube.CubeList()
 mean_SH_Ex_10yr=iris.cube.CubeList()
 stdev_Tr_10yr=iris.cube.CubeList()
 stdev_NH_Ex_10yr=iris.cube.CubeList()
 stdev_SH_Ex_10yr=iris.cube.CubeList()

 for i in range (10, 200+1, 10):
   st=i-10
   mean_Tr_10yr.append( (tropics_MO[st:i] - tropics_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_Tr_10yr.append( (tropics_MO[st:i] - tropics_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
   mean_NH_Ex_10yr.append( (NH_Ex_MO[st:i] - NH_Ex_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_NH_Ex_10yr.append( (NH_Ex_MO[st:i] - NH_Ex_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
   mean_SH_Ex_10yr.append( (SH_Ex_MO[st:i] - SH_Ex_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_SH_Ex_10yr.append( (SH_Ex_MO[st:i] - SH_Ex_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
  #print mean_Gl_10yr
 #print mean_Gl_10yr[0]
 #print mean_Gl_10yr[-1]

 ##30 years differences 
 mean_Tr_30yr=iris.cube.CubeList()
 mean_NH_Ex_30yr=iris.cube.CubeList()
 mean_SH_Ex_30yr=iris.cube.CubeList()
 stdev_Tr_30yr=iris.cube.CubeList()
 stdev_NH_Ex_30yr=iris.cube.CubeList()
 stdev_SH_Ex_30yr=iris.cube.CubeList()

 for i in range (30, 200+1, 30):
   st=i-30
   mean_Tr_30yr.append( (tropics_MO[st:i] - tropics_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_Tr_30yr.append( (tropics_MO[st:i] - tropics_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
   mean_NH_Ex_30yr.append( (NH_Ex_MO[st:i] - NH_Ex_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_NH_Ex_30yr.append( (NH_Ex_MO[st:i] - NH_Ex_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
   mean_SH_Ex_30yr.append( (SH_Ex_MO[st:i] - SH_Ex_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_SH_Ex_30yr.append( (SH_Ex_MO[st:i] - SH_Ex_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )

 ##50 years differences 
 mean_Tr_50yr=iris.cube.CubeList()
 mean_NH_Ex_50yr=iris.cube.CubeList()
 mean_SH_Ex_50yr=iris.cube.CubeList()
 stdev_Tr_50yr=iris.cube.CubeList()
 stdev_NH_Ex_50yr=iris.cube.CubeList()
 stdev_SH_Ex_50yr=iris.cube.CubeList()

 for i in range (50, 200+1, 50):
   st=i-50
   mean_Tr_50yr.append( (tropics_MO[st:i] - tropics_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_Tr_50yr.append( (tropics_MO[st:i] - tropics_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
   mean_NH_Ex_50yr.append( (NH_Ex_MO[st:i] - NH_Ex_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_NH_Ex_50yr.append( (NH_Ex_MO[st:i] - NH_Ex_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
   mean_SH_Ex_50yr.append( (SH_Ex_MO[st:i] - SH_Ex_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_SH_Ex_50yr.append( (SH_Ex_MO[st:i] - SH_Ex_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )


 ##100 years differences 
 mean_Tr_100yr=iris.cube.CubeList()
 mean_NH_Ex_100yr=iris.cube.CubeList()
 mean_SH_Ex_100yr=iris.cube.CubeList()
 stdev_Tr_100yr=iris.cube.CubeList()
 stdev_NH_Ex_100yr=iris.cube.CubeList()
 stdev_SH_Ex_100yr=iris.cube.CubeList()

 for i in range (100, 200+1, 100):
   st=i-100
   mean_Tr_100yr.append( (tropics_MO[st:i] - tropics_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_Tr_100yr.append( (tropics_MO[st:i] - tropics_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
   mean_NH_Ex_100yr.append( (NH_Ex_MO[st:i] - NH_Ex_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_NH_Ex_100yr.append( (NH_Ex_MO[st:i] - NH_Ex_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )
   mean_SH_Ex_100yr.append( (SH_Ex_MO[st:i] - SH_Ex_AR[st:i]).collapsed('time', iris.analysis.MEAN) )
   stdev_SH_Ex_100yr.append( (SH_Ex_MO[st:i] - SH_Ex_AR[st:i]).collapsed('time', iris.analysis.STD_DEV) )


 ##200 years differences 
 mean_Tr_200yr= (tropics_MO - tropics_AR).collapsed('time', iris.analysis.MEAN) 
 stdev_Tr_200yr= (tropics_MO - tropics_AR).collapsed('time', iris.analysis.STD_DEV) 
 mean_NH_Ex_200yr=  (NH_Ex_MO - NH_Ex_AR).collapsed('time', iris.analysis.MEAN) 
 stdev_NH_Ex_200yr= (NH_Ex_MO - NH_Ex_AR).collapsed('time', iris.analysis.STD_DEV) 
 mean_SH_Ex_200yr=  (SH_Ex_MO - SH_Ex_AR).collapsed('time', iris.analysis.MEAN) 
 stdev_SH_Ex_200yr= (SH_Ex_MO - SH_Ex_AR).collapsed('time', iris.analysis.STD_DEV) 

 mean_tropics=iris.cube.CubeList([mean_Tr_10yr.merge_cube(), mean_Tr_30yr.merge_cube(), mean_Tr_50yr.merge_cube(), mean_Tr_100yr.merge_cube(), mean_Tr_200yr])
 stdev_tropics=iris.cube.CubeList([stdev_Tr_10yr.merge_cube(), stdev_Tr_30yr.merge_cube(), stdev_Tr_50yr.merge_cube(), stdev_Tr_100yr.merge_cube(), stdev_Tr_200yr])
 mean_NH_Ex=iris.cube.CubeList([mean_NH_Ex_10yr.merge_cube(), mean_NH_Ex_30yr.merge_cube(), mean_NH_Ex_50yr.merge_cube(), mean_NH_Ex_100yr.merge_cube(), mean_NH_Ex_200yr])
 stdev_NH_Ex=iris.cube.CubeList([stdev_NH_Ex_10yr.merge_cube(), stdev_NH_Ex_30yr.merge_cube(), stdev_NH_Ex_50yr.merge_cube(), stdev_NH_Ex_100yr.merge_cube(), stdev_NH_Ex_200yr])
 mean_SH_Ex=iris.cube.CubeList([mean_SH_Ex_10yr.merge_cube(), mean_SH_Ex_30yr.merge_cube(), mean_SH_Ex_50yr.merge_cube(), mean_SH_Ex_100yr.merge_cube(), mean_SH_Ex_200yr])
 stdev_SH_Ex=iris.cube.CubeList([stdev_SH_Ex_10yr.merge_cube(), stdev_SH_Ex_30yr.merge_cube(), stdev_SH_Ex_50yr.merge_cube(), stdev_SH_Ex_100yr.merge_cube(), stdev_SH_Ex_200yr])
 '''
 #save in a csv file mean and stdev
 fileout=open(fname, "w") #'a'
 fileout.write("mean_10yr_tropics" + "stdev_10yr_tropics")
 fileout.write(str(mean_tropics[0].data) + str(stdev_tropics[0].data))
 fileout.write("mean_10yr_NH_Ex" + "stdev_10yr_NH_Ex")
 fileout.write(str(mean_NH_Ex[0].data) + str(stdev_NH_Ex[0].data))
 fileout.write("mean_10yr_SH_Ex" + "stdev_10yr_SH_Ex")
 fileout.write(str(mean_SH_Ex[0].data) + str(stdev_SH_Ex[0].data))
 fileout.write("mean_30yr_tropics" + "stdev_30yr_tropics")
 fileout.write(str(mean_tropics[1].data) + str(stdev_tropics[1].data))
 fileout.write("mean_30yr_NH_Ex" + "stdev_30yr_NH_Ex")
 fileout.write(str(mean_NH_Ex[1].data) + str(stdev_NH_Ex[1].data))
 fileout.write("mean_30yr_SH_Ex" + "stdev_30yr_SH_Ex")
 fileout.write(str(mean_SH_Ex[1].data) + str(stdev_SH_Ex[1].data))
 fileout.write("mean_50yr_tropics" + "stdev_50yr_tropics")
 fileout.write(str(mean_tropics[2].data) + str(stdev_tropics[2].data))
 fileout.write("mean_50yr_NH_Ex" + "stdev_50yr_NH_Ex")
 fileout.write(str(mean_NH_Ex[2].data) + str(stdev_NH_Ex[2].data))
 fileout.write("mean_50yr_SH_Ex" + "stdev_50yr_SH_Ex")
 fileout.write(str(mean_SH_Ex[2].data) + str(stdev_SH_Ex[2].data))
 fileout.write("mean_100yr_tropics" + "stdev_100yr_tropics")
 fileout.write(str(mean_tropics[3].data) + str(stdev_tropics[3].data))
 fileout.write("mean_100yr_NH_Ex" + "stdev_100yr_NH_Ex")
 fileout.write(str(mean_NH_Ex[3].data) + str(stdev_NH_Ex[3].data))
 fileout.write("mean_100yr_SH_Ex" + "stdev_100yr_SH_Ex")
 fileout.write(str(mean_SH_Ex[3].data) + str(stdev_SH_Ex[3].data))
 fileout.write("mean_200yr_tropics" + "stdev_200yr_tropics")
 fileout.write(str(mean_tropics[4].data) + str(stdev_tropics[4].data))
 fileout.write("mean_200yr_NH_Ex" + "stdev_200yr_NH_Ex")
 fileout.write(str(mean_NH_Ex[4].data) + str(stdev_NH_Ex[4].data))
 fileout.write("mean_200yr_SH_Ex" + "stdev_200yr_SH_Ex")
 fileout.write(str(mean_SH_Ex[4].data) + str(stdev_SH_Ex[4].data))
 '''

 return mean_tropics, stdev_tropics, mean_NH_Ex, stdev_NH_Ex, mean_SH_Ex, stdev_SH_Ex

#######################################################################################################################################################################################

########################################### Plotting function ##########################################################################################################################

def Plot(MO, tropics_MO, tropics_AR, NH_Ex_MO, NH_Ex_AR, SH_Ex_MO, SH_Ex_AR, mean_tropics, mean_NH_Ex, mean_SH_Ex, title1="title", title2="title", title3="title", title4="title", ylab1="lab", ylab2="lab"):

 ##create time coordinate for plots: 
 years=MO.shape[0]
 #times= pd.date_range(start='1850-01-01', periods=years, freq='AS') 
 times=np.arange(0,years,1)

 #plotting 
 fig=plt.figure(figsize=(20,10))

 #subplot1
 ax = fig.add_subplot(221) 
 ax=plt.gca()
 ax.set_title('Tropics')
 #fig.autofmt_xdate()

 plt.plot( times, tropics_MO.data, c='grey',  linewidth=2, label='MONSOON' )
 plt.plot( times, tropics_AR.data, c='black', linestyle= '--', linewidth=2, label='ARCHER')
 plt.title( title1, fontsize=16)
 plt.xlabel('Time (years)', fontsize=14)
 plt.ylabel( ylab1, fontsize=14)
 plt.ticklabel_format(axis='y',useOffset=False)
 #ax.xaxis.set_major_locator(mdates.YearLocator(20))
 ax.add_artist(AnchoredText('a', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )
 plt.legend(loc=2)



 ##subplot2
 ax = fig.add_subplot(222) 
 ax=plt.gca() 
 ax.set_title('NH Ex')
 #fig.autofmt_xdate()

 plt.plot( times, NH_Ex_MO.data, c='grey',  linewidth=2, label='MONSOON' )
 plt.plot( times, NH_Ex_AR.data, c='black', linestyle= '--', linewidth=2, label='ARCHER')
 plt.title( title2, fontsize=16)
 plt.xlabel('Time (years)', fontsize=14)
 plt.ylabel( ylab1, fontsize=14)
 plt.ticklabel_format(axis='y',useOffset=False)
 #ax.xaxis.set_major_locator(mdates.YearLocator(20))
 ax.add_artist(AnchoredText('b', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )
 plt.legend(loc=2)


 ##subplot3
 ax = fig.add_subplot(223)
 ax=plt.gca()
 ax.set_title('SH Ex')
 #fig.autofmt_xdate()

 plt.plot( times, SH_Ex_MO.data, c='grey',  linewidth=2, label='MONSOON' )
 plt.plot( times, SH_Ex_AR.data, c='black', linestyle= '--', linewidth=2, label='ARCHER')
 plt.title( title3 , fontsize=16)
 plt.xlabel('Time (years)', fontsize=14)
 plt.ylabel( ylab1 , fontsize=14)
 #plt.ticklabel_format(useOffset=False)
 plt.ticklabel_format(axis='y',useOffset=False)
 #ax.xaxis.set_major_locator(mdates.YearLocator(20))
 ax.add_artist(AnchoredText('c', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) )
 plt.legend(loc=2)


 ##subplot4
 ax = fig.add_subplot(224)
 ax=plt.gca()
 ax.set_title('Delta')
 ax.axhline(y=0 , color='black', linestyle='--') #plot zero line

 #t_scales=np.array([ [10]*20, [30]*6, [50]*4, [100]*2, [200] ])

 for n in range (0, 4): 
  for i in range(0, mean_NH_Ex[n].shape[0]):
   nh=plt.scatter(n-0.15, mean_NH_Ex[n][i].data, marker='D',  s=60, facecolors='none', edgecolors='orange' )
 plt.scatter(4-0.15, mean_NH_Ex[4].data, marker='D',  s=60, facecolors='none', edgecolors='orange'  )
 
 for n in range (0, 4): 
  for i in range(0, mean_tropics[n].shape[0]):
   gl=plt.scatter(n, mean_tropics[n][i].data, marker='o', s=60, facecolors='none', edgecolors='green')
 plt.scatter(4, mean_tropics[4].data, marker='o', s=60, facecolors='none', edgecolors='green')


 for n in range (0, 4): 
  for i in range(0,mean_SH_Ex[n].shape[0]):
   sh=plt.scatter(n+0.15, mean_SH_Ex[n][i].data, c='red' , marker='x',  s=60)
 plt.scatter(4+0.15, mean_SH_Ex[4].data, c='red', marker='x',  s=60 )


 plt.title( title4 , fontsize=16)
 plt.xlabel('Timescale (years)', fontsize=14)
 plt.ylabel(ylab2, fontsize=14)
 plt.legend((gl, nh, sh), ('Tropics', 'NH Extratropics', 'SH Extratropics'), loc=1, fontsize=13)
 plt.xticks([0, 1, 2, 3, 4], [ 10,  30,  50,  100,  200])
 ax.add_artist(AnchoredText('d', loc=2, borderpad=0.0, prop=dict(fontsize=16), bbox_to_anchor=(1.,1.), bbox_transform=ax.transAxes ) ) 
 
# ##################################################################################################################################################################################

#load SW data 
MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_TOA_LW&SW_18502050.nc',  'toa_outgoing_shortwave_flux' ) 
AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_TOA_LW&SW_18502050.nc',     'toa_outgoing_shortwave_flux')
#Now use G_N_S function on SW data:
SW_tropics_MO, SW_tropics_AR, SW_NH_Ex_MO, SW_NH_Ex_AR, SW_SH_Ex_MO, SW_SH_Ex_AR = G_N_S(MO, AR)
#Now use MO_minus_AR function on SW data to compute differences, associated stats and save it in a file:
SW_mean_tropics, SW_stdev_tropics, SW_mean_NH_Ex, SW_stdev_NH_Ex, SW_mean_SH_Ex, SW_stdev_SH_Ex = MO_minus_AR(SW_tropics_MO, SW_tropics_AR, SW_NH_Ex_MO, SW_NH_Ex_AR, SW_SH_Ex_MO, SW_SH_Ex_AR, "stats_SW.csv")
#call the plotting function:
Plot(MO, SW_tropics_MO, SW_tropics_AR, SW_NH_Ex_MO, SW_NH_Ex_AR, SW_SH_Ex_MO, SW_SH_Ex_AR, SW_mean_tropics, SW_mean_NH_Ex, SW_mean_SH_Ex, title1="TOA outgoing SW flux in the tropics", title2="TOA outgoing SW flux in the NH extratropics", title3="TOA outgoing SW flux in the SH extratropics", title4="SW TOA differences as a function of the averaging timescale",  ylab1="SW TOA (W/m$^{2}$)", ylab2="SW$_{MO}$-SW$_{AR}$ (W/m$^{2}$)")


#load LW data 
MO=iris.load_cube('/nerc/n02/n02/vittoria/uar766_TOA_LW&SW_18502050.nc',  'toa_outgoing_longwave_flux' ) 
AR=iris.load_cube('/nerc/n02/n02/vittoria/uas245_TOA_LW&SW_18502050.nc',     'toa_outgoing_longwave_flux')
#Now use G_N_S function on LW data:
LW_tropics_MO, LW_tropics_AR, LW_NH_Ex_MO, LW_NH_Ex_AR, LW_SH_Ex_MO, LW_SH_Ex_AR = G_N_S(MO, AR)
#Now use MO_minus_AR function on LW data to compute differences, associated stats and save it in a file:
LW_mean_tropics, LW_stdev_tropics, LW_mean_NH_Ex, LW_stdev_NH_Ex, LW_mean_SH_Ex, LW_stdev_SH_Ex = MO_minus_AR(LW_tropics_MO, LW_tropics_AR, LW_NH_Ex_MO, LW_NH_Ex_AR, LW_SH_Ex_MO, LW_SH_Ex_AR, "stats_LW.csv")
#call the plotting function:
Plot(MO, LW_tropics_MO, LW_tropics_AR, LW_NH_Ex_MO, LW_NH_Ex_AR, LW_SH_Ex_MO, LW_SH_Ex_AR, LW_mean_tropics, LW_mean_NH_Ex, LW_mean_SH_Ex, title1="TOA outgoing LW flux in the tropics", title2="TOA outgoing LW flux in the NH extratropics", title3="TOA outgoing LW flux in the SH extratropics", title4="LW TOA differences as a function of the averaging timescale",  ylab1="LW TOA (W/m$^{2}$)", ylab2="LW$_{MO}$-LW$_{AR}$ (W/m$^{2}$)")


plt.show()


