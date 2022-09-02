# OGIVE_finder
script for identifying OGIVES from digital elevation data

### OGIvE_finder.py

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

### OGIVE_functions.py 

This script contains various gee function to convert linestrings to usable points in order to compute zonal stats.
