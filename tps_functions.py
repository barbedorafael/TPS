# =============================================================================
# Functions
# =============================================================================

# =============================================================================
# Compute stream mask
# =============================================================================

from osgeo import gdal, osr
import numpy as np

def flow_ds(rowL, colL, fdr):
    
    dir1 = fdr[rowL,colL] == 1
    dir2 = fdr[rowL,colL] == 2
    dir3 = fdr[rowL,colL] == 4
    dir4 = fdr[rowL,colL] == 8
    dir5 = fdr[rowL,colL] == 16
    dir6 = fdr[rowL,colL] == 32
    dir7 = fdr[rowL,colL] == 64
    dir8 = fdr[rowL,colL] == 128
    
    rowL[dir1] = rowL[dir1]
    colL[dir1] = colL[dir1] + 1
    
    rowL[dir2] = rowL[dir2] + 1
    colL[dir2] = colL[dir2] + 1
    
    rowL[dir3] = rowL[dir3] + 1
    colL[dir3] = colL[dir3]
    
    rowL[dir4] = rowL[dir4] + 1
    colL[dir4] = colL[dir4] - 1
    
    rowL[dir5] = rowL[dir5]
    colL[dir5] = colL[dir5] - 1
    
    rowL[dir6] = rowL[dir6] - 1
    colL[dir6] = colL[dir6] - 1

    rowL[dir7] = rowL[dir7] - 1
    colL[dir7] = colL[dir7]
    
    rowL[dir8] = rowL[dir8] - 1
    colL[dir8] = colL[dir8] + 1
    
    return rowL, colL

def _pop_rim(data, nodata=0):
    data[:,0] = nodata
    data[:,-1] = nodata
    data[0,:] = nodata
    data[-1,:] = nodata
    return data

def get_str_mask(fdr, fac, tpi, thmax, a, th_tpi):
    
    fdr = _pop_rim(fdr)
    fac_th = np.where(tpi<=th_tpi, a, thmax)

    # Stream definition
    mask = np.where((fac >= fac_th) | (np.isin(fdr, [0, -1, 255])), 1, 0)
    
    # Follow fdir to fill "cutted" values
    strIdx = np.where((mask == 1) & (fac <= thmax))
    flow_ds(strIdx[0], strIdx[1], fdr)
    noSDS = mask[strIdx]==0
    strIdx = (strIdx[0][noSDS],   # Takes only arguments where downstream has no stream
               strIdx[1][noSDS])
    
    while np.any(mask[strIdx] == 0):
        mask[strIdx] = 1
        flow_ds(strIdx[0], strIdx[1], fdr)


    return mask

# Save raster
def array2tif(outname, raster, mx):
    # writing output 
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(outname, mx.shape[1], mx.shape[0], 1, gdal.GDT_Float32)
    ds.SetProjection(raster.GetProjection())
    ds.SetGeoTransform(raster.GetGeoTransform())
    ds.GetRasterBand(1).WriteArray(mx)
    ds = None

# =============================================================================
# Convert raster to vector line ! to be optimized
# =============================================================================

from shapely import geometry as geo
import shapely.vectorized
import fiona
from fiona.crs import from_epsg

def pixelOffset2coord(raster,xOffset,yOffset):
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    coordX = originX+pixelWidth*xOffset+pixelWidth*0.5
    coordY = originY+pixelHeight*yOffset-pixelWidth*0.5
    return coordX, coordY

def str_nbhd(strl, fdr, X, Y):
    dir1 = Y, X+1    # right
    dir2 = Y+1, X+1  # down-right
    dir3 = Y+1, X    # down
    dir4 = Y+1, X-1  # down-left
    dir5 = Y, X-1    # left
    dir6 = Y-1, X-1  # up-left
    dir7 = Y-1, X    # up
    dir8 = Y-1, X+1  # up-right
    
    # assign neighborhood value of fdr map
    # assign neighborhood value of stream map
    # assign value of fdr that would drain to cell (upstream->downstream)
    dir_dict = {
                1: (fdr[dir1], strl[dir1], 16),
                2: (fdr[dir2], strl[dir2], 32),
                3: (fdr[dir3], strl[dir3], 64),
                4: (fdr[dir4], strl[dir4], 128),
                5: (fdr[dir5], strl[dir5], 1),
                6: (fdr[dir6], strl[dir6], 2),
                7: (fdr[dir7], strl[dir7], 4),
                8: (fdr[dir8], strl[dir8], 8)
                }
    
    # headwater = 0, middle-stream = 1, junction = 2+
    strup = 0
    for dr in dir_dict.values():
        if dr[1] == 1 and dr[0]==dr[2]:
            strup += 1
    
    return strup

def arrayClip(rasref, array, poly):
    
    src = rasref
    ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
    lrx = ulx + (src.RasterXSize * xres)
    lry = uly + (src.RasterYSize * yres)
    
    x  = np.linspace(ulx, lrx, src.RasterXSize)
    y  = np.linspace(uly, lry, src.RasterYSize)
    
    xx,yy = np.meshgrid(x,y)
    
    inpoly = shapely.vectorized.contains(poly,  xx,yy)    
    
    if array.dtype == np.dtype('float32'):
        array[~inpoly] = np.nan
    else:
        array[~inpoly] = 0
    
    return array
    
def array2shp(rasterFdr, mask, fdr, fac, lyr_fn = "strLayer.shp", pixelValue=1):
    
    # array2dict
    count = 0
    strList = np.where(mask == pixelValue)
    nbhdList = np.zeros((strList[0].size), dtype=np.uint8)
    for indexY in strList[0]:
        indexX = strList[1][count]
        
        if np.any((indexY<10) | (indexY>fdr.shape[0]-10) |
                  (indexX<10) | (indexX>fdr.shape[1]-10)):
            nbhdList[count] = 2
        else:
            nbhdList[count] = (str_nbhd(mask, fdr, indexX, indexY))
        count += 1
    
    # Loop to write MultiLineString
    headList = np.where(nbhdList==0)[0]
    stretchIn = []
    stretches = []
    for i in headList:
        Y = strList[0][i]
        X = strList[1][i]
        Xcoord, Ycoord = pixelOffset2coord(rasterFdr,X,Y)
        point = [Xcoord, Ycoord]
        path = 1
        while path == 1:
            point0 = point
            line = [point]
            idx = [(Y,X)]
            acc = [fac[Y,X]]
            stretch = 1
            while stretch == 1:
                Y, X = flow_ds(np.array(Y), np.array(X), fdr)
                Xcoord, Ycoord = pixelOffset2coord(rasterFdr,X,Y)
                point = [Xcoord, Ycoord]         
                
                line.append(point)
                idx.append((Y,X))
                acc.append(fac[Y,X])
                nbhdVal = nbhdList[(strList[0]==Y) & (strList[1]==X)]
                

                if np.any((Y<10) | (Y>fdr.shape[0]-10) |
                          (X<10) | (X>fdr.shape[1]-10) |
                          (fdr[Y,X] == 0) | (fdr[Y,X] == 255) |
                          (not nbhdVal)):
                    stretch = 0
                    path = 0                    

                elif nbhdVal > 1:
                    stretch = 0
            
            if fdr[Y,X] == 255:
                continue
            else:
                stretches.append(geo.LineString(line))
                stretchIn.append(point0)          
            if point in stretchIn:
                path = 0
    
    # Save to file
    schema = {
        'geometry': 'LineString',
        'properties': {'id': 'int'},
    }
    
    with fiona.open(lyr_fn, "w", driver="ESRI Shapefile",
                    crs=from_epsg(4326), schema=schema) as layer:
        for i, stretch in enumerate(stretches):
            layer.write({
                'geometry': {
                    'type': 'LineString',
                    'coordinates': stretch
                },
                'properties': {'id': i},
            })
    
    return


