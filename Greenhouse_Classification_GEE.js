Map.addLayer(geometry) 

/**
* Function to mask clouds using the Sentinel-2 QA band
* @param {ee.Image} image Sentinel-2 image
* @return {ee.Image} cloud masked Sentinel-2 image
*/
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 3;
  var cirrusBitMask = 1 << 5;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000);
}


// Map the function over one year of data and take the median.
// Load Sentinel-2 TOA reflectance data.
var image_1 = ee.ImageCollection('COPERNICUS/S2_SR')
                  .filterBounds(geometry)
                  .filterDate('2019-03-01', '2019-05-31')
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                  .map(maskS2clouds)
                  .select('B8','B4','B2','B3','B11').median();
var image_2 = ee.ImageCollection('COPERNICUS/S2_SR')
                  .filterBounds(geometry)
                  .filterDate('2019-06-01', '2019-08-31')
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                  .map(maskS2clouds)
                  .select('B8','B4','B2','B3','B11').median();
var image_3 = ee.ImageCollection('COPERNICUS/S2_SR')
                  .filterBounds(geometry)
                  .filterDate('2019-09-01', '2019-10-31')
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                  .map(maskS2clouds)
                  .select('B8','B4','B2','B3','B11').median();
var image_4 = ee.ImageCollection('COPERNICUS/S2_SR')
                  .filterBounds(geometry)
                  .filterDate('2019-01-01', '2019-02-28')
                  // Pre-filter to get less cloudy granules.
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                  .map(maskS2clouds)
                  .select('B8','B4','B2','B3','B11').median();

var rgbVis = {
  min: 0.0,
  max: 0.3,
  gamma: 1.4,
  bands: ['B4', 'B3', 'B2'],
};
Map.addLayer(image_1, rgbVis, 'RGB_0305');
Map.addLayer(image_2, rgbVis, 'RGB_0608');
Map.addLayer(image_3, rgbVis, 'RGB_0910');
Map.addLayer(image_4, rgbVis, 'RGB_1201');


var compall = ee.Image.cat([image_1,image_3]);
var newimage0 = compall.reproject({crs: compall.projection().crs(), scale: 10});
var sentinel2 = newimage0.clip(geometry);

//DEM dataset                  
var DEMdataset = ee.Image('JAXA/ALOS/AW3D30_V1_1').clip(geometry);
var elevation = DEMdataset.select('AVE');
var elevationVis = {
  min: 0.0,
  max: 4000.0,
  palette: ['0000ff', '00ffff', 'ffff00', 'ff0000', 'ffffff'],
};


//slope dataset
var slope = ee.Terrain.slope(elevation).clip(geometry);

//GLCM
var getL8GLCM = function(data) {
  var asm = data.int32().glcmTexture({size: 4}).select('B2_asm');
  var contrast = data.int32().glcmTexture({size: 4}).select('B2_contrast');
  var corr = data.int32().glcmTexture({size: 4}).select('B2_corr');
  var vari = data.int32().glcmTexture({size: 4}).select('B2_var');
  var idm = data.int32().glcmTexture({size: 4}).select('B2_idm');
  var savg = data.int32().glcmTexture({size: 4}).select('B2_savg');
  var ent = data.int32().glcmTexture({size: 4}).select('B2_ent');
  var diss = data.int32().glcmTexture({size: 4}).select('B2_diss');
  return data.addBands(asm).addBands(contrast).addBands(corr).addBands(vari).addBands(idm)
            .addBands(ent).addBands(savg).addBands(diss);
};

//Spectral Index
function cal_ndvi(data){
  var ndvi0 = data.expression('(nir-red)/(nir+red)',{'nir' :data.select('B8'),'red' :data.select('B4')});
  var ndvi = ee.Image(ndvi0);
  return ndvi;
}
function cal_ndbi(data){
  var ndbi0 = data.expression('(swir1-nir)/(nir+swir1)',{'nir' :data.select('B8'),'swir1' : data.select('B11')});
  var ndbi = ee.Image(ndbi0);
  return ndbi;
}
function cal_mndwi(data){
  var ndwi0 = data.expression('(green-mir)/(mir+green)',{'mir' :data.select('B11'),'green' :data.select('B3')});
  var mndwi = ee.Image(ndwi0);
  return mndwi;
}
function cal_savi(data){
  var savi0 = data.expression('1.5*(nir-red)/(nir+red+0.5)',{'nir' :data.select('B8'),'red' :data.select('B4')});
  var savi = ee.Image(savi0);
  return savi;
}

var ndvi1 = cal_ndvi(image_1).rename('ndvi1');
var ndvi2 = cal_ndvi(image_2).rename('ndvi2');
var ndvi3 = cal_ndvi(image_3).rename('ndvi3');
var ndvi4 = cal_ndvi(image_4).rename('ndvi4');

var ndbi1 = cal_ndbi(image_1).rename('ndbi1');
var ndbi2 = cal_ndbi(image_2).rename('ndbi2');
var ndbi3 = cal_ndbi(image_3).rename('ndbi3');
var ndbi4 = cal_ndbi(image_4).rename('ndbi4');

var mndwi1 = cal_mndwi(image_1).rename('mndwi1');
var mndwi2 = cal_mndwi(image_2).rename('mndwi2');
var mndwi3 = cal_mndwi(image_3).rename('mndwi3');
var mndwi4 = cal_mndwi(image_4).rename('mndwi4');

var savi1 = cal_savi(image_1).rename('savi1');
var savi2 = cal_savi(image_2).rename('savi2');
var savi3 = cal_savi(image_3).rename('savi3');
var savi4 = cal_savi(image_4).rename('savi4');

var glcm_1 = getL8GLCM(image_1).select('B2_asm','B2_contrast','B2_corr','B2_var','B2_idm','B2_savg','B2_ent','B2_diss');
var glcm_2 = getL8GLCM(image_2).select('B2_asm','B2_contrast','B2_corr','B2_var','B2_idm','B2_savg','B2_ent','B2_diss');
var glcm_3 = getL8GLCM(image_3).select('B2_asm','B2_contrast','B2_corr','B2_var','B2_idm','B2_savg','B2_ent','B2_diss');
var glcm_4 = getL8GLCM(image_4).select('B2_asm','B2_contrast','B2_corr','B2_var','B2_idm','B2_savg','B2_ent','B2_diss');

Map.centerObject(image_1, 9);
var ndviParams = {min: -1, max: 1, palette: ['blue', 'white', 'green']};
Map.addLayer(ndvi1, ndviParams, 'NDVI1 image_1');
Map.addLayer(ndbi1, ndviParams, 'NDBI1 image_1');

Map.addLayer(ndvi2, ndviParams, 'NDVI2 image_2');
Map.addLayer(ndbi2, ndviParams, 'NDBI2 image_2');

Map.addLayer(ndvi3, ndviParams, 'NDVI3 image_3');
Map.addLayer(ndbi3, ndviParams, 'NDBI3 image_3');

Map.addLayer(ndvi4, ndviParams, 'NDVI4 image_4');
Map.addLayer(ndbi4, ndviParams, 'NDBI4 image_4');

// Create 1000 random points in the region.
var randomPointsDP = ee.FeatureCollection.randomPoints(DP,5000,0);
var randomPointsNONDP = ee.FeatureCollection.randomPoints(NONDP,5000,0);
var addproperty0 = function(feature) {
  return feature.set('landcover',0);
};
var addproperty1 = function(feature) {
  return feature.set('landcover',1);
};

// Map the area getting function over the FeatureCollection.
var pointnondp  = randomPointsNONDP.map(addproperty0);
var pointdp = randomPointsDP.map(addproperty1);


var pointnondp = pointnondp.randomColumn();
var pointdp = pointdp.randomColumn();

// Split Training data and Validation data
var split = 0.7;  // 7:3
// DP
var Training_data_dp = pointdp.filter(ee.Filter.lte('random', split));
print('Training_data_dp',Training_data_dp)
var Validation_data_dp = pointdp.filter(ee.Filter.gte('random', split));
print('Validation_data_dp',Validation_data_dp)
// NoDP
var Training_data_nodp = pointnondp.filter(ee.Filter.lte('random', split));
print('Training_data_nodp',Training_data_nodp)
var Validation_data_nodp = pointnondp.filter(ee.Filter.gte('random', split));
print('Validation_data_nodp',Validation_data_nodp)
//merging points
var Training_points = Training_data_dp.merge(Training_data_nodp);
var Validation_points = Validation_data_dp.merge(Validation_data_nodp);
// build dataset for Train
var dataset = image_1.addBands(image_3).addBands(glcm_1).addBands(glcm_3)
                      .addBands(ndvi1).addBands(ndbi1).addBands(mndwi1).addBands(savi1)
                      .addBands(ndvi3).addBands(ndbi3).addBands(mndwi3).addBands(savi3)
                      .addBands(elevation).addBands(slope);
print('DATASAT:',dataset);
// Sample the input imagery to get a FeatureCollection of training data.
var training = dataset.sampleRegions({
  collection: Training_points,
  properties: ['landcover'], 
  scale:30
});
// Make a Random Forest classifier and train it.
var classifier = ee.Classifier.smileRandomForest(150,4).train(training,'landcover');

// Classify the image.
var classified = dataset.classify(classifier);
var kernel = ee.Kernel.square({radius: 1});
// Perform an erosion followed by a dilation, display.
var opened = classified
            .focal_min({kernel: kernel, iterations: 1})
            .focal_max({kernel: kernel, iterations: 1});

// Define a palette for the DP classification.
var palette = [
'#FFFFFF',//  nonDP (0) // 
'#F08080', // DP (1)  //                                    

];
Map.setCenter(119.0079, 37.0765, 9);
Map.addLayer(classified, {min: 0, max: 1, palette: palette}, 'DP Classification');
Map.addLayer(opened, {min: 0, max: 1, palette: palette}, 'DP OPENED');

//Sample the input imagery to get a FeatureCollection of Validation data.
var Validation_points = dataset.sampleRegions({
  collection: Validation_points,
  properties: ['landcover'],
  scale: 30
});

// Classify the validation data.
var validated = Validation_points.classify(classifier,'testlandcover');
 
// Get a confusion matrix representing expected accuracy.
var testAccuracy = validated.errorMatrix('landcover','testlandcover');
print('Validation error matrix: ', testAccuracy);
print('Validation overall accuracy: ', testAccuracy.accuracy());

//Export image and table 
var exportAccuracy1 = ee.Feature(null, {matrix: testAccuracy.array()});
var exportAccuracy2 = ee.Feature(null, {accuracy:testAccuracy.accuracy()});



Export.table.toDrive({
  collection: ee.FeatureCollection(exportAccuracy2),
  description: '2019_DP_acc_WeiF_T20210610',
  fileFormat: 'CSV'
});

Export.image.toDrive({
  image:classified,
  description: "2019_CALSS_WeiF_T20210610",
  fileNamePrefix: "2019_CALSS_WeiF_T20210610",
  scale: 30,
  region:geometry,
  maxPixels: 1e13
});
// Export.image.toDrive({
//   image:opened,
//   description: "2019_OPENED_WeiF_T20210610",
//   fileNamePrefix: "2019_OPENED_WeiF_T20210610",
//   scale: 30,   //when scale = 10 . error user memory limit exceeded
//   region:geometry,
//   maxPixels: 1e13
// });

