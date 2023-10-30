def centroid_extract(feature):

  """ function for extracting centroids """

  return feature.centroid(**{'maxError': 1}).select([]).set('land_cover', 1)

def mask_clouds(image):

  """ function for removing clouds """

  qa = image.select('QA60')
  cloudBitMask = 1 << 10
  cirrusBitMask = 1 << 11

  mask = qa.bitwiseAnd(cloudBitMask).eq(0).And(qa.bitwiseAnd(cirrusBitMask).eq(0))

  return image.updateMask(mask).divide(10000)


# Function for calculating vegetation indices
def vi_calculaton(image):

  """ function for calculating vegetation indices """

  NDVI = image.normalizedDifference(['B8', 'B4']).rename('NDVI') #Normalized difference vegetation index

    #Normalized difference water index
  NDWI = image.expression(
    '(GREEN - NIR)/(GREEN + NIR)', {
     'NIR': image.select('B8'),
     'GREEN': image.select('B3')
    }).rename('NDWI')

    #Normalized difference bare index
  NDBI = image.expression(
    '(SWIR - NIR)/(SWIR + NIR)', {
     'SWIR': image.select('B11'),
     'NIR': image.select('B8')
    }).rename('NDBI')

    #Green index (G.I.)
  GI = image.expression(
    '(NIR / Green)', {
     'Green': image.select('B3'),
     'NIR': image.select('B8')
    }).rename('GI')


    #Chlorophyll index (CI)
  CI = image.expression(
    '(NIR / RedEdge) - 1', {
     'NIR': image.select('B8'),
     'RedEdge': image.select('B5')
    }).rename('CI')

    #Green chlorophyll vegetation index (GCVI)
  GCVI = image.expression(
    'NIR / Green - 1', {
     'Green': image.select('B3'),
     'NIR': image.select('B8')
    }).rename('GCVI')

    #Land surface water index (LSWI)
  LSWI = image.expression(
    '(NIR - SWIR)/(NIR + SWIR)', {
     'SWIR': image.select('B11'),
     'NIR': image.select('B8')
    }).rename('LSWI')

    #Enhanced vegetation index (EVI)
  EVI = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('B8'),
      'RED': image.select('B4'),
      'BLUE': image.select('B2')
        }).toFloat().rename('EVI')

    # Bare soil index (BSI)
  BSI = image.expression(
    '(SWIR + Red) - (NIR + Blue) / (SWIR + Red) + (NIR + Blue)', {
      'NIR': image.select('B8'),
      'SWIR': image.select('B11'),
      'Red': image.select('B4'),
      'Blue': image.select('B2')
        }).toFloat().rename('BSI')

  return image.addBands([NDVI, NDWI, NDBI, GI, GCVI, LSWI, EVI, BSI, CI])


def normalize(image):

  """ function for normalizing images """

  bandNames = image.bandNames()

    #/ Compute min and max of the image
  minDict =  image.reduceRegion(**{
      'reducer' : ee.Reducer.min(),
      'geometry' : AOI_1,
      'scale' : 10,
      "maxPixels" : 1e9,
      'bestEffort' : True,
      "tileScale" : 16
    })

  maxDict =  image.reduceRegion(**{
    'reducer' : ee.Reducer.max(),
    'geometry' : AOI_1,
    'scale' : 10,
    'maxPixels' : 1e9,
    'bestEffort' : True,
    'tileScale' : 16
    })

  mins = ee.Image.constant(minDict.values(bandNames))
  maxs = ee.Image.constant(maxDict.values(bandNames))

  normalized =  image.subtract(mins).divide(maxs.subtract(mins))

  return normalized


def monthly_composite(start_date, end_date):
    
    #Variables to retain when returning the layertacked image
    selected_feature_names = ['et','temp','LSWI','vhIwAscMean','NDWI','GCVI','NDBI','elevation','NDVI']
    
    # cropland mask from Dynamic World
    global_lc2020 = ee.ImageCollection('ESA/WorldCover/v200') # With Dynamic World, you can get monthly landcover maps
    #Create a mask of cropland (from ESA)
    global_cropland_lc2020 = global_lc2020.first().eq(40).clipToCollection(AOI_1).selfMask() # Extract the cropland areas from the ESA land cover. This will be used as mask

    # Variables to consider: ['et','temp','LSWI','vhIwAscMean','NDWI','GCVI','NDBI','elevation','NDVI']
    ######################################################### Sentinel-2 ####################################
    S2img = ee.ImageCollection("COPERNICUS/S2_SR") \
        .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',10)).filterBounds(AOI_1).filterDate(start_date,      end_date).map(mask_clouds).median().select('B.*')  # updateMsk instead of clip

    # Calculate Vegetation indices
    #NDVI, NDWI, NDBI,G.I.,GCVI,LSWI, EVI, BSI, Terrain elevation,Terrain slope (still more to be added)
    
    S2_vi = vi_calculaton(S2img)
    S2_vi = S2_vi.select(['LSWI','NDWI','GCVI','NDBI','NDVI'])
    
    ######################################################### End of Sentinel-2 ####################################

    ######################################################### Sentinel-1 ####################################
    S1img = ee.ImageCollection('COPERNICUS/S1_GRD') \
                          .filterDate(start_date, end_date).filterBounds(AOI_1)
    
    #Filter the Sentinel-1 collection by metadata properties.
    vvVhIw = S1img.filter(ee.Filter.listContains('transmitterReceiverPolarisation',   'VV')).filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')).filter(ee.Filter.eq('instrumentMode', 'IW'));
    
    #Separate ascending and descending orbit images into distinct collections.
    vvVhIwAsc = vvVhIw.filter(
      ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));
    vvVhIwDesc = vvVhIw.filter(
      ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));
    
    # Calculate temporal means for various observations to use for visualization.
    # Mean VH ascending.
    vhIwAscMean = vvVhIwAsc.select('VH').mean();
    # Mean VH descending.
    vhIwDescMean = vvVhIwDesc.select('VH').mean();
    # Mean VV for combined ascending and descending image collections.
    vvIwAscDescMean = vvVhIwAsc.merge(vvVhIwDesc).select('VV').mean();
    #Mean VH for combined ascending and descending image collections.
    vhIwAscDescMean = vvVhIwAsc.merge(vvVhIwDesc).select('VH').mean();
    
    S1Img = ee.ImageCollection.fromImages([vhIwAscMean]).toBands();
    
    ######################################################### End of Sentinel-1 ####################################
    
    ######################################################### Climatic data ####################################
    # prec = ee.ImageCollection('UCSB-CHG/CHIRPS/DAILY') \
                             # .filter(ee.Filter.date(start_date, end_date)).select('precipitation').mean()  #precipitation product: Chirp
    
    et = ee.ImageCollection("FAO/WAPOR/2/L1_AETI_D") \
                             .filter(ee.Filter.date(start_date, end_date)).mean()  ## evapotranspiration product: Wapor
    
    temp = ee.ImageCollection('MODIS/061/MOD11A1').select('LST_Day_1km').filter(ee.Filter.date(start_date, end_date)).mean()  ## temperature product:MOD11A1.061 Terra Land Surface Temperature and Emissivity Daily Global 1km
    

    clim = ee.ImageCollection.fromImages([et, temp]).toBands()
    ######################################################### End of Climatic data ####################################
    
    ######################################################### Terrain variables ####################################
    # I'm using the copernicus DEM
    glo30 = ee.ImageCollection("projects/sat-io/open-datasets/GLO-30");
    
    dem = glo30.mosaic().setDefaultProjection('EPSG:3857',None,30)
    
    #Elevation
    dem_clip = dem.clip(AOI_1)  #Clip the dem to AOI
    
    #Slope
    # slope = ee.Terrain.slope(dem_clip);
    
    dem_slope = ee.ImageCollection.fromImages([dem_clip]).toBands()
    
    
    ######################################################### End of Terrain variables ####################################
    ######################################################### Combine all the variables ####################################
    all_variables = ee.ImageCollection.fromImages([S2_vi, S1Img, dem_slope, clim]).toBands().updateMask(global_cropland_lc2020)
    
    img_for_classification = all_variables.select([
                                               '0_LSWI','0_NDWI','0_GCVI','0_NDBI','0_NDVI','1_0_VH','2_0_b1','3_0_L1_AETI_D','3_1_LST_Day_1km']).rename('LSWI', 'NDWI','GCVI', 'NDBI',  
                                           'NDVI','vhIwAscMean','elevation','et','temp'
                                        )   #rename bands in the composite image
    
    final_img = img_for_classification.select(selected_feature_names)
    return final_img
    