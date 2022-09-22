# -*- coding: utf-8 -*-
"""
Created on Thursday sep 15th 10:41:24 2022

@author: Yoram

OGIVE_finder


This script pulls the elevation and slope profiles of ArcticDEM along a centerline given as a shapefile 
by the user. Use earth engine to easily draw one. 
Arcitc DEM is pulled using Google earth engine tools.

The pulled elevation profile is detrended by a smoothed version of the same profile. 

Then, we run a fast fourrier transform to get at the periodicity of the data. This seems to work best on the elevation vs. the slope profile.

This script also pulls the average ice velocity from the Millan estimations, present as RGI_v1-6 in earth engine assets (not great)

Running a periodogram rather than an fft allows to retain information on the location along the centerline of the periodicity, 
which allows to identify the location of the OGIVES.
Ogives are flagged when the amplitude of oscillation at the wavelength corresponding to ~1 year of displacement (avg. ice velocities from Millan et
                                                                                                                 et al. 2021)
is more than 1.5 m on the DEM.  
"""


# import modules and functions
import numpy as np 
from functools import reduce
import ee as ee
import geopandas as gp 
import geemap 
import os
from os.path import exists
from numpy import genfromtxt
from scipy.fft import fft, fftfreq
from scipy.signal import spectrogram

import matplotlib.pyplot as plt


# OGIVE_ functions contains the gee script to get at points 
from OGIVE_functions import line_to_points, buffer_points

# for GEE
ee.Initialize()

#%% work on slope profile or elev profile ?
work_on_slope = False

#%% read in data 
print('Loading data...')

#%% select a glacier #####################################################################
RGIglacier = ee.FeatureCollection('GLIMS/current').filter('glac_name=="Kennicott Glacier"')

# optional to use ID for non named glaciers 
 # RGIglacier = ee.FeatureCollection('GLIMS/current').filter('glac_id=="G209330E63184N"')


# glacier name for file out
glacier_name = 'Gilkey'

########################################################################################


#%% load ArcticDEM, and convert to slope (in degrees)
elev = ee.Image("UMN/PGC/ArcticDEM/V3/2m_mosaic")#.select('elevation') 
slope = ee.Terrain.slope(elev)

#%% load in ice thickness from the GEE asset library 
# this could be improved to make the script more flexible...)
RGI = ee.Image('projects/ee-yoramterleth/assets/RGI_AK3')


#%% calculate the bed elevation and slope  
bed_elev = elev.subtract(RGI) 
bed_slope = ee.Terrain.slope(bed_elev)


#%% centerlines from Zhang 2022

#read in the shapefile 
#centerlines = gp.read_file('C:/Users/Yoram/OneDrive - University of Idaho\Desktop/PhD pos/SURGE_CYCLES/peclet_nb/python_gee/centerlines/01_rgi60_Alaska_final_lines.shp')
centerlines = gp.read_file('C:/Users/Yoram/OneDrive - University of Idaho\Desktop/PhD pos/OGIVES/findr/centerlines/cl_gates.shp')
#centerlines = gp.read_file('C:/Users/Yoram/OneDrive - University of Idaho\Desktop/PhD pos/OGIVES/findr/centerlines/RGI/AK_centerlines.shp')

# set the CRS to coincide with the ARctic dem one
centerlines = centerlines.to_crs(crs='EPSG:4269')

#%% select the centerlines to a specific glacier 

# get the outline from earth engine 
GLIMS = RGIglacier# ee.FeatureCollection('GLIMS/current').filter('glac_name=="Haenke Glacier"'); 

# convert this to a shapefile that can be used with geopandas
geemap.ee_to_shp(GLIMS,'C:/Users/Yoram/OneDrive - University of Idaho\Desktop/PhD pos/SURGE_CYCLES/peclet_nb/python_gee/centerlines/'+glacier_name+'_outline.shp' )

GLIMS_shp = gp.read_file('C:/Users/Yoram/OneDrive - University of Idaho\Desktop/PhD pos/SURGE_CYCLES/peclet_nb/python_gee/centerlines/'+glacier_name+'_outline.shp' )
GLIMS_shp = GLIMS_shp.to_crs(crs='EPSG:4269')

cl = gp.clip(centerlines, GLIMS_shp)

#%% plot the glacier outline and the centerlines we are pulling
ax = GLIMS_shp.plot()
ax = cl.plot()

#%% load ArcticDEM, and convert to slope (in degrees)
elev = ee.Image("UMN/PGC/ArcticDEM/V3/2m_mosaic") 
slope = ee.Terrain.slope(elev)

#%% get the strings from the centerlines and convert to gee features 
print('Converting centerlines to earth engine features...')


# import shapely

# # Put the sub-line coordinates into a list of sublists
# outcoords = [list(i.coords) for i in cl.geometry[1457]]

# # Flatten the list of sublists and use it to make a new line
# outline = shapely.geometry.LineString([i for sublist in outcoords for i in sublist])


#%% 

centerline = cl.geometry[0] # [1113] #[1282] #outline #merged_line #cl.geometry[15653][0] #[36][0]
x,y = centerline.coords.xy 
cords = np.dstack((x,y)).tolist()
double_list = reduce(lambda x,y: x+y, cords)
ee_centerline = ee.Geometry.LineString(double_list)

#% transorm the gee linestring to points 
print('Getting points from centerlines...')

# nb of points we want to pull
point_nb = 1000 

# length of the ee string
length = ee_centerline.length()

# step interval at which to get points
step = ee_centerline.length().divide(point_nb)

# distances at which each point needs to appear
distances = ee.List.sequence(0,length, step)

# cut up the linestring into pieces of that length 
lines= ee.Geometry(ee_centerline).cutLines(distances).geometries()  

# clear the pointz variable just in case 
pointz = []   

# apply the line to points function  
pointz = line_to_points(ee.Geometry(ee_centerline), 1000)

# apply the bufferpoints function 
pts_topo = pointz #.map(buffer_points)

# make a buffered points collection too: we will use this pulled data to detrend our profile later
pts_buffered = pointz.map(buffer_points)


#% pull the actual data 

# concatenate the variables of interest from the dem as bands of an image, and add a dummy system start time 
topo = ee.Image.cat(elev, slope, RGI).set('system:time_start', ee.Date('2000-01-01').millis())
topo = topo.rename(['elev', 'slope', 'vel'])


# set up the csv step needed for geemap.zonal_statistics()... this is clunky but it works remarkably well

# give the filenames
outdir = 'C:/Users/Yoram/OneDrive - University of Idaho\Desktop/PhD pos/OGIVES/findr/centerlines/'
outfile = outdir +'zonalstats4.csv'
bufferfile = outdir + 'bufferstats4.csv'

# remove the files if they already exists
if exists(outfile):
    os.remove(outfile)
    
if exists(bufferfile):
    os.remove(bufferfile)

# run zonal statistics functions, and save profiles to csv

print('Running zonal stats...')

geemap.zonal_statistics(topo,pts_topo,outfile,statistics_type='MEAN',scale=2)

geemap.zonal_statistics(topo,pts_buffered,bufferfile,statistics_type='MEAN',scale=2)



#%finding periodicity in the elevation profiles:
print('Finding periodicity in profile...')

# get the non buffered profile 
data = genfromtxt(outfile, delimiter=',')

# pull it apart into components 
zz = data[1:,0] 
xx = data[1:,5]
yy = data[1:,4]
ss = data[1:,1]
offset = data[1:,6]
ind = data[1:,3]
vel = data[1:,2]


# get smoothed profile and pull the relevant parts out
bufferdata = genfromtxt(bufferfile, delimiter=',')
zz_smooth = bufferdata[1:,0]
ss_smooth = bufferdata[1:,1] 
offset_smooth = bufferdata[1:,3]


#% detrend the profiles 
zz_small_scale = zz-zz_smooth
ss_small_scale = ss-ss_smooth


#% select elev or slope to work on 
if work_on_slope:
    targ_data = ss_small_scale
    targ_true = ss
else: 
    targ_data = zz_small_scale 
    targ_true = zz
    

#% COMPUTE FFT

# (this second part can be tested with a centerline from code editor for which we know there are ogives, before automating. )

# remmove nans 
zf = [x for x in targ_data if np.isnan(x) == False]

# romeove nans from centerline data
offset_f = offset[np.isnan(targ_data) == False]



# number of sample points 
N = len(zf)

# sample spacing 
T = np.mean(np.diff(offset))

# run fft 
yf = fft(zf)

# get frequencies
xf = fftfreq(N, T)[0:N//2]
xw = 1/xf

#% COMPUTE SPECTROGRAM
f, t, Sxx = spectrogram(np.array(targ_data),1/T,window=('tukey', 0.25),nperseg=100,noverlap=70,scaling='spectrum')

# convert output freuqnecies to wavelength
wavelength = 1/f 

#%% VISUALIZE
print('Plotting...')

# plot the fft
fig1, ax1 = plt.subplots(figsize= (10,5))
plt.plot(xw, 2.0/N * np.abs(yf[0:N//2]))
plt.grid()
ax1.set_ylabel('fft')
ax1.set_xlabel('wavelenght (m)')
plt.show()



# plot the spectrogram 
fig3, ax3 = plt.subplots(figsize=(10,5))
aa= plt.pcolormesh(t, wavelength[1:], Sxx[1:])#,vmin=0,vmax=200) #, shading='flat')
plt.ylabel('Oscillation Wavelength (m)')
plt.xlabel('Distance along centerline (m)')

ax3.set_ylim(20,500)

    

plt.show()


#%% get the velocity of interest depending on location along centerline
target_wavelength = vel /1 

select_array = np.zeros(len(t))

window_margin = 3 # in steps 

# select from Sxx at the target_wavelength, -+ 50 m 
for i in range(len(t)):
    
    # find index of the target wavelenght
    targ_loc = np.abs(t[i] - offset).argmin()
    ind_cent = (np.abs(wavelength[1:] - target_wavelength[targ_loc])).argmin()
    
    # select a slightly wider area when possible 
    indices = np.arange(ind_cent - window_margin, ind_cent + window_margin, 1).tolist()
    
    # select the Sxx value in that area
    try:
        select_array[i] = np.amax(Sxx[indices,i])
    except:
        if indices[-1] > len(Sxx[:,0]):
            select_array[i] = np.amax(Sxx[indices[0]:,i])
            print('Reduced window at upper array bound.At d=' + str(t[i]))
        elif indices[0] < 0:
            select_array[i] = np.amax(Sxx[:indices[-1],i])
            print('Reduced window at lower array bound. At d=' + str(t[i]))
        else: 
            print('Error: indexing should not fail. Skipping.At d=' + str(t[i]))
            continue

    
#%% select the coordinates for the ogives 

# need to implementn this, interp the max osci amplitude 
# at the spatial res of the offset vecotrs, and make the trheshold ynamic...
interp_select_array = np.interp(offset, t, select_array)

# select area that is above 1 m, or above percentile threshold...
select_xx = xx[interp_select_array>1.5] #np.percentile(interp_select_array,85)]
select_yy = yy[interp_select_array>1.5] #np.percentile(interp_select_array,85)]

print('plotting')


# plot the detended profile used
fig2, (ax2a,ax2,ax2b,ax2c) = plt.subplots(4,1,figsize=(12,13))
ax2.plot(offset,targ_data)
ax2c.set_xlabel('distance along centerline (m)')

ax2a.plot(offset,targ_true)

ax2b.plot(offset, vel)
ax2b.set_ylabel('Average ice velocity (m/year)')

ax2c.plot(offset,interp_select_array)
ax2c.set_ylabel('Max. osc. Amp. - 1 year period (m)')

if work_on_slope:
    ax2.set_ylabel('detrended slope (deg)')
    fig3.colorbar(aa,label='Oscillation amplitude (slope)')
    ax2a.set_ylabel('slope (deg)')
else:
    ax2.set_ylabel('detrended elev. m')
    fig3.colorbar(aa,label='Oscillation amplitude (m)')
    ax2a.set_ylabel('elev. (m a.s.l.)')
    

fig3, ax3 = plt.subplots(figsize=(10,10))
GLIMS_shp.plot(ax=ax3,alpha=0.7)
ax3.scatter(select_xx,select_yy, c='red')
cl.plot(ax=ax3,color='green',linewidth=2) 
ax3.set_xlabel('lon (deg)')
ax3.set_ylabel('lat (deg)')

print('Done!')



