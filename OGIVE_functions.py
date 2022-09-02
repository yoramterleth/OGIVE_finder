# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:31:08 2022

@author: Yoram

ogive helper functions

"""

#%% declare packages
# import numpy as np 
# from functools import reduce
import ee as ee
# import geopandas as gp 
# import pandas as pd

# import matplotlib.pyplot as plt

# from functools import partial 
#from list_all_assets import get_asset_list

#from ee_converter import feature2ee
#from ee_converter import feature2ee

ee.Initialize()



#%% make point feature, actually makes the points

def makePointFeature(coord, offset):
    pt = ee.Algorithms.GeometryConstructors.Point(coord);
    coords = pt.coordinates(); # set coords as feature property
    
    newpt = ee.Feature(pt).set({'offset':offset,
                                'lat':ee.List(coords).get(0),
                                'lon':ee.List(coords).get(1)})
    return newpt


#%% point line function, determines where the points go

def point_line(s):
    line = ee.List(s).get(0)
    offset = ee.List(s).get(1)
    # make first point 
    point = makePointFeature(ee.Geometry(line).coordinates().get(0), offset)
    return point


#%% line to points function, converts a slinetring into individual sampling points

def line_to_points(linestring, count):
    length = linestring.length()
    step = linestring.length().divide(count)
    distances = ee.List.sequence(0,length,step)
    
    lines = ee.Geometry(linestring).cutLines(distances).geometries()
    
    points = lines.zip(distances).map(point_line)
      
    pts = ee.FeatureCollection(points)
    
    return pts



#%% buffer points function, makes feature collection that contains plygons around each point of area given here
# will be mapped over all points in collection 

def buffer_points(pt):
    radius = 50 
    bp = ee.Feature(pt).buffer(radius)
    
    return bp
