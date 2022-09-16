// Mini script to draw in a quick centerline and export it as a shapefile. Also helpful to figure out 
// specific glacier names or IDs in the RGI, using the AOI option 
////////////////////////////////////////////////////////////////////////////////////////////

// import the glims outlines
var GLIMS = ee.FeatureCollection('GLIMS/current').filter('glac_name=="Kennicott Glacier"'); 

// uncoment to look for glaciers within user drawn AOI:
// var GLIMS = ee.FeatureCollection('GLIMS/current').filter(ee.Filter.bounds(AOI)) ; 

print(GLIMS)
Map.addLayer(GLIMS,{color: 'purple',alpha: 0.5},'Glacier centerline')

// Export the FeatureCollection to a KML file.
Export.table.toDrive({
  collection: ee.FeatureCollection(cl_gates),
  description:'cl_gates',
  fileFormat: 'SHP'
});
