//hyperparameter
var start_1 = '2019-01-01';
var end_1 = '2019-04-01';
var start_2 = '2019-06-01';
var end_2 = '2019-10-01';
var rgbVis = {min: 0.0, max: 0.3, gamma: 1.4, bands: ['B4', 'B3', 'B2'],};
var ndviParams = {min: -1, max: 1, palette: ['blue', 'white', 'green']};
var elevationVis = {min: 0.0, max: 4000.0, palette: ['0000ff', '00ffff', 'ffff00', 'ff0000', 'ffffff']};

// --------------------------------------------------------------------------------------------------------------------------------
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

// Function to load remote sensing images
// Map the function over one year of data and take the median.
// Load Sentinel-2 TOA reflectance data.
function load_data(start, end, table_geo){
  // var table_geo = table_geo;
  var image = ee.ImageCollection('COPERNICUS/S2_SR')
                    .filterBounds(table_geo)
                    .filterDate(start, end)
                    // Pre-filter to get less cloudy granules.
                    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                    .map(maskS2clouds)
                    .select('B1','B2','B3','B4','B8','B11','B12').median();
                    
  return image;
  
}

var image_1 =  load_data(start_1, end_1, table);   
var image_2 =  load_data(start_2, end_2, table); 
Map.addLayer(image_1.clip(table), rgbVis, 'beijing_spring--4_1');
Map.addLayer(image_2.clip(table), rgbVis, 'beijing_spring--4_2');

// Building Feature Function
// --spectral index
function cal_ndvi(data){
  var ndvi0 = data.expression('(nir-red)/(nir+red)',{'nir':data.select('B8'),'red':data.select('B4')});
  var ndvi = ee.Image(ndvi0);
  return ndvi;
}
function cal_ndbi(data){
  var ndbi0 = data.expression('(swir1-nir)/(nir+swir1)',{'nir':data.select('B8'),'swir1': data.select('B11')});
  var ndbi = ee.Image(ndbi0);
  return ndbi;
}
function cal_mndwi(data){
  var ndwi0 = data.expression('(green-swir1)/(swir1+green)',{'swir1':data.select('B11'),'green':data.select('B3')});
  var mndwi = ee.Image(ndwi0);
  return mndwi;
}
function cal_savi(data){
  var savi0 = data.expression('1.5*(nir-red)/(nir+red+0.5)',{'nir':data.select('B8'),'red':data.select('B4')});
  var savi = ee.Image(savi0);
  return savi;
}
function cal_apgi(data){
  var apgi0 = data.expression('100*aerosols*red*(2*nir-red-swir2)/(2*nir+red+swir2)',
              {'aerosols':data.select('B1'),'red':data.select('B4'),'nir':data.select('B8'),'swir2':data.select('B12')});
  var apgi = ee.Image(apgi0);
  return apgi;
}
function cal_pghi(data){
  var pghi0 = data.expression('blue/swir2',{'blue': data.select('B2'),'swir2':data.select('B12')});
  var pghi = ee.Image(pghi0);
  return pghi;
}
function cal_pgi(data){
  var pgi0 = data.expression('(100*(blue*(nir-red)))/(1-(blue+green+nir))/3',
  {'blue': data.select('B2'),'nir':data.select('B8'),'red':data.select('B4'),'green':data.select('B3')});
  var pgi = ee.Image(pgi0);
  return pgi;
}
function cal_pmli(data){
  var pmli0 = data.expression('(swir1-red)/(swir1+red)',
  {'swir1':data.select('B11'), 'red':data.select('B4')});
  var pmli = ee.Image(pmli0);
  return pmli;
}
function cal_rpgi(data){
  var rpgi0 = data.expression('100*blue/(1-(nir+blue+green))/3',
  {'blue':data.select('B2'), 'nir':data.select('B8'), 'green':data.select('B3')});
  var rpgi = ee.Image(rpgi0);
  return rpgi;
}

// --GLCM
function get_GLCM(data){
  var asm = data.int32().glcmTexture({size: 4}).select('B2_asm');
  var contrast = data.int32().glcmTexture({size: 4}).select('B2_contrast');
  var corr = data.int32().glcmTexture({size: 4}).select('B2_corr');
  var vari = data.int32().glcmTexture({size: 4}).select('B2_var');
  var idm = data.int32().glcmTexture({size: 4}).select('B2_idm');
  var savg = data.int32().glcmTexture({size: 4}).select('B2_savg');
  var ent = data.int32().glcmTexture({size: 4}).select('B2_ent');
  var diss = data.int32().glcmTexture({size: 4}).select('B2_diss');
  
  return asm.addBands(contrast).addBands(corr).addBands(vari).addBands(idm)
            .addBands(ent).addBands(savg).addBands(diss);
}

// Calculating Feature map
var ndvi1 = cal_ndvi(image_1).rename('ndvi1');
var ndbi1 = cal_ndbi(image_1).rename('ndbi1');
var mndwi1 = cal_mndwi(image_1).rename('mndwi1');
var savi1 = cal_savi(image_1).rename('savi1');
var apgi1 = cal_apgi(image_1).rename('apgi1');
var pghi1 = cal_pghi(image_1).rename('pghi1');
var pgi1 = cal_pgi(image_1).rename('pgi1');
var pmli1 = cal_pmli(image_1).rename('pmli1');
var rpgi1 = cal_rpgi(image_1).rename('rpgi1');

var ndvi2 = cal_ndvi(image_2).rename('ndvi2');
var ndbi2 = cal_ndbi(image_2).rename('ndbi2');
var mndwi2 = cal_mndwi(image_2).rename('mndwi2');
var savi2 = cal_savi(image_2).rename('savi2');
var apgi2 = cal_apgi(image_2).rename('apgi2');
var pghi2 = cal_pghi(image_2).rename('pghi2');
var pgi2 = cal_pgi(image_2).rename('pgi2');
var pmli2 = cal_pmli(image_2).rename('pmli2');
var rpgi2 = cal_rpgi(image_2).rename('rpgi2');

var glcm_1 = get_GLCM(image_1).select('B2_asm','B2_contrast','B2_corr','B2_var','B2_idm','B2_savg','B2_ent','B2_diss');
var glcm_2 = get_GLCM(image_2).select('B2_asm','B2_contrast','B2_corr','B2_var','B2_idm','B2_savg','B2_ent','B2_diss');

// add property to each class
var NonDP_a = ee.FeatureCollection(NonDP);
var DP_a = ee.FeatureCollection(DP);

var addproperty0 = function(feature) {return feature.set('landcover',0);};
var addproperty1 = function(feature) {return feature.set('landcover',1);};

var Point_NonDP_b = NonDP_a.map(addproperty0); 
var Point_DP_b = DP_a.map(addproperty1); 

// split training data and validation data
var NonDP_c = Point_NonDP_b.randomColumn();
var DP_c = Point_DP_b.randomColumn();

var split = 0.7;
var training_NonDP_c = NonDP_c.filter(ee.Filter.lt('random', split));            
var validation_NonDP_c = NonDP_c.filter(ee.Filter.gte('random', split));    

var training_DP_c = DP_c.filter(ee.Filter.lt('random', split));            
var validation_DP_c = DP_c.filter(ee.Filter.gte('random', split)); 

var training_mix = ee.FeatureCollection(training_NonDP_c).merge(training_DP_c);
var validation_mix = ee.FeatureCollection(validation_NonDP_c).merge(validation_DP_c);

print('traing_mix',training_mix.size());
print('validation_mix',validation_mix.size());

// Filter out the null property values and try again.
var training = training_mix.filter(
  ee.Filter.notNull(training_mix.first().propertyNames())
);
var validation = validation_mix.filter(
  ee.Filter.notNull(validation_mix.first().propertyNames())
);

// Building Feature Dataset for training operation
var dataset = image_1.addBands(image_2)
              .addBands(glcm_1).addBands(glcm_2)
              .addBands(ndvi1).addBands(ndbi1).addBands(mndwi1).addBands(savi1).addBands(apgi1).addBands(pghi1).addBands(pgi1).addBands(pmli1).addBands(rpgi1)
              .addBands(ndvi2).addBands(ndbi2).addBands(mndwi2).addBands(savi2).addBands(apgi2).addBands(pghi2).addBands(pgi2).addBands(pmli2).addBands(rpgi2);
print("dataset",dataset);              

var dataset_clip = dataset.clip(table);
print("dataset_clip",dataset_clip);


// var geometry = ee.FeatureCollection('users/mrsuperniu/WeiFang_GEE_Grid_40km').filter(ee.Filter.eq('Id',grid_id)).geometry();
var geometry = table;
Map.addLayer(geometry); 
// Overlay the points on the imagery to get training.
// get each spectral index feature value with training data point 
var training_RF = dataset_clip.sampleRegions({
  collection: training,
  properties: ['landcover'], 
  scale:10,
  tileScale: 8 
});
print('training_RF:',training_RF);

var classifier = ee.Classifier.smileRandomForest(150,4).train(training_RF,'landcover');
// Classify the image.
var classified = dataset_clip.clip(geometry).classify(classifier); // 'clip' here is necessary
// var kernel = ee.Kernel.square({radius: 1});
// Perform an erosion followed by a dilation, display.
// var opened = classified
//           .focal_min({kernel: kernel, iterations: 1})
//           .focal_max({kernel: kernel, iterations: 1});

var palette = [
'#ffffff',    // 'NonDP'
'#00f5ff',    // 'DP'
];
Map.addLayer(classified, {min: 0, max: 1, palette: palette}, 'DP Classification');
// Map.addLayer(opened.clip(table), {min: 0, max: 1, palette: palette}, 'DP OPENED');

var validation_RF = dataset_clip.sampleRegions({
  collection: validation,
  properties: ['landcover'],
  scale: 10,
  tileScale: 8
});
print('validation_RF:',validation_RF);
// Classify the validation data.
var validated = validation_RF.classify(classifier,'testlandcover');
// Get a confusion matrix representing expected accuracy.
var testAccuracy = validated.errorMatrix('landcover','testlandcover');
print('Validation error matrix: ', testAccuracy);
print('Validation overall accuracy: ', testAccuracy.accuracy());
