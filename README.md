# Topographic Position based Stream definition (TPS)
Python code for stream definition in a Digital Elevation Model (DEM) conditioned by Topographic Position Index (TPI).

#### Reference:
Rafael Barbedo, Vinícius Siqueira & Walter Collischonn (2022). [Topographic Position-based Stream definition (TPS): A simple method to address spatial variability of drainage density in stream networks. Hydrological Sciences Journal.](https://www.tandfonline.com/doi/full/10.1080/02626667.2022.2047190)

### Inputs required
- Digital Elevation Model (DEM)
- Flow direction (FDR) - must be in [HydroSHEDS](https://developers.google.com/earth-engine/datasets/catalog/WWF_HydroSHEDS_03DIR#:~:text=HydroSHEDS%20is%20based%20on%20elevation,vary%20from%201%20to%20128.) format (1=E, 2=SE, 4=S, 8=SW, 16=W, 32=NW, 64=N, 128=NE)
- Flow accumulation (FAC) - in Area OR cells
- Topographic Position Index (TPI) - from [Weiss (2001)](http://www.jennessent.com/downloads/tpi-poster-tnc_18x22.pdf); a radius of 900m works well for this method with a 90m resolution DEM (~10 DEM cells)

All these maps are easily obtained using a GIS processing tool in a DEM (*e.g.* QGIS).
You can also download them via Google Earth Engine. Just [click here](https://code.earthengine.google.com/f38773265937b3f731020af260679492) to check the GEE code (using MERIT Hydro database).

### Parameters
- TPI threshold (V)
- Min Area threshold (Amin)
- Max Area threshold (Amax)

Stream is defined where `FAC > Amin AND TPI < V` OR `FAC > Amax`. Then, gaps are filled using FDR.

![image](https://user-images.githubusercontent.com/83959435/119422091-98f2bc00-bcd6-11eb-8bc3-9e919b3767f3.png | width=100)


## Example usage
1. Run the `tps_functions.py` file in your python IDLE

2. Load the input rasters
```python
# Source file to copy saving format
region = 'Corrente'
srcfn = 'GIS/'+region+'_dem_merit_90.tif'

# Load rasters
dem = gdal.Open(srcfn).ReadAsArray()
fdr = gdal.Open('GIS/'+region+'_fdr_merit_90.tif').ReadAsArray()
fac = gdal.Open('GIS/'+region+'_fac_merit_90.tif').ReadAsArray()
tpi = gdal.Open('GIS/'+region+'_tpi_merit_90.tif').ReadAsArray()
```
3. Define parameters
```python
# V (m)
v = -8 

# Amin (km²)
amin = 1.6

# Amax (km²)
amax = 320
```

4. Compute stream network with TPS method
```python
str_mask = get_str_mask(fdr, fac, tpi, amax, amin, v)

# Save tif file
array2tif('GIS/'+region+'_str_mask.tif', srcfn, str_mask)
```

5. \(Optional) Convert stream mask to line vector
```python
str_line = array2gdf(rasref, mask_tpi, fdr, fac)

# Save shp file
str_line.to_file('GIS/'+region+'_str_line.shp')
```


