# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 10:41:24 2022

@author: Yoram

OGIVE_finder


This script pulls the elevation and slope profiles of ArcticDEM along a centerlines from 
Zhang et al (2022) - A new global dataset of mountain glacier centerlines and lengths (https://doi.org/10.5194/essd-14-3889-2022); 
available here: https://doi.org/10.11922/sciencedb.01643 

This dataset needs to be availaible in the working directory.
Arcitc DEM is pulled using Google earth engine tools.

The pulled elevation profile is detrended by a smoothed version of the same profile. 

Then, we run a fast fourrier transform to get at the periodicity of the data. This seems to work best on the slope profile.

Ultimetely, we should be able to compare with ice velocity from Itslive, and ensure that the periodicity is annual. 
Running a periodogram rather than an fft should also allow us to retain information on the location along the centerline of the periodicity, 
which should allow us to identify the location of the OGIVES. 
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

import matplotlib.pyplot as plt


# OGIVE_ functions contains the gee script to get at points 
from OGIVE_functions import line_to_points, buffer_points

# for GEE
ee.Initialize()

#%% work on slope profile or elev profile ?
work_on_slope = True 

#%% read in data 
print('Loading data...')

# read in the shapefile 
centerlines = gp.read_file('C:/Users/Yoram/OneDrive - University of Idaho\Desktop/PhD pos/OGIVES/findr/centerlines/01_rgi60_Alaska_final_lines.shp')



# select only the really long glaciers for now 
centerlines = centerlines[centerlines.MaxL>30000]

# set the CRS to coincide with the ARctic dem one
centerlines = centerlines.to_crs(crs='EPSG:4269')

# now get a specific centerline for testing
centerline = centerlines.geometry[12][0]

#%% load ArcticDEM, and convert to slope (in degrees)
elev = ee.Image("UMN/PGC/ArcticDEM/V3/2m_mosaic") 
slope = ee.Terrain.slope(elev)

#%% get the strings from the centerlines and convert to gee features 
print('Converting centerlines to ear engine features...')

x,y = centerline.coords.xy 
cords = np.dstack((x,y)).tolist()
double_list = reduce(lambda x,y: x+y, cords)
ee_centerline = ee.Geometry.LineString(double_list)

#%% transorm the gee linestring to points 
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


#%% pull the actual data 

# concatenate the variables of interest from the dem as bands of an image, and add a dummy system start time 
topo = ee.Image.cat(elev, slope).set('system:time_start', ee.Date('2000-01-01').millis())


# set up the csv step needed for geemap.zonal_statistics()... this is clunky but it works remarkably well

# give the filenames
outdir = 'C:/Users/Yoram/OneDrive - University of Idaho\Desktop/PhD pos/OGIVES/findr/centerlines/'
outfile = outdir +'zonalstats3.csv'
bufferfile = outdir + 'bufferstats.csv'

# remove the files if they already exists
if exists(outfile):
    os.remove(outfile)
    
if exists(bufferfile):
    os.remove(bufferfile)

# run zonal statistics functions, and save profiles to csv

print('Running zonal stats...')

geemap.zonal_statistics(topo,pts_topo,outfile,statistics_type='MEAN',scale=2)

geemap.zonal_statistics(topo,pts_buffered,bufferfile,statistics_type='MEAN',scale=2)



#%% finding periodicity in the elevation profiles:
print('Finding periodicity in profile...')

# get the non buffered profile 
data = genfromtxt(outfile, delimiter=',')

# pull it apart into components 
zz = data[1:,0] 
xx = data[1:,3]
yy = data[1:,4]
ss = data[1:,1]
offset = data[1:,2]
ind = data[1:,5]

# get smoothed profile and pull the relevant parts out
bufferdata = genfromtxt(bufferfile, delimiter=',')
zz_smooth = bufferdata[1:,0]
ss_smooth = bufferdata[1:,1] 
offset_smooth = bufferdata[1:,2]


#%% detrend the profiles 
zz_small_scale = zz-zz_smooth
ss_small_scale = ss-ss_smooth


#%%  select elev or slope to work on 
if work_on_slope:
    data = ss_small_scale
else: 
    data = zz_small_scale 
    

#%% COMPUTE FFT

# (this second part can be tested with a centerline from code editor for which we know there are ogives, before automating. )

# remmove nans 
zf = [x for x in data if np.isnan(x) == False]


# number of sample points 
N = len(zf)

# sample spacing 
T = np.mean(np.diff(offset))

# run fft 
yf = fft(zf)

# get frequencies
xf = fftfreq(N, T)[0:N//2]

#%% VISUALIZE
print('Plotting...')

# plot the fft
fig1, ax1 = plt.subplots(figsize= (10,5))
plt.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
plt.grid()
ax1.set_ylabel('fft')
ax1.set_xlabel('period')
plt.show()

# plot the detended profile used
fig2, ax2 = plt.subplots(figsize=(10,5))
plt.plot(zf)
ax2.set_ylabel('detrended slope')
ax2.set_xlabel('distance along centerline (m)')


print('Done!')



