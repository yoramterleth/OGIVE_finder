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

Ultimately, we should be able to compare with ice velocity from Itslive, and ensure that the periodicity is annual. 
Running a periodogram rather than an fft should also allow us to retain information on the location along the centerline of the periodicity, 
which should allow us to identify the location of the OGIVES. 

### OGIVE_findr_single_glacier_v2.py

This script does the same as above, but on a specific user defined glacier. Abandoned the Zhang profiles in favor of those that come with the RGI database, as they are much more continuous and still cover a lot of tributaries. Can input a specific shapefile centerline for testing. Now also pulls average spatially distributed velocity from Millan et al 2021 estimations, further steps are to hone in on the ogives in the spectrogram, and to implement it automatically over a lot of centerlines. 

### OGIVE_findr_single_glacier_v2.py

This script does the same as above, but it loops over all the centerlines that are in the RGI database for a specific glacier. The glacier selection happens where the user inputs the name (or id) of the GLIMS glacier. 

### OGIVE_functions.py 

This script contains various gee function to convert linestrings to usable points in order to compute zonal stats.
