#### Just want to clarify again here, `okada85.m` is a creation from [okada85.m](https://github.com/IPGP/deformation-lib/tree/master/okada). My Matlab codes are designed to leverage its functionality.  

---
All Matlab variables are the format of ***structure***. It should be rather straight-forward to modify datasets by replacing data with the correct field.  
4 main data structures are used and created throughout the slip inversion:  
* Displacement structure (e.g. DataStruct)
* Fault model structure (e.g. FaultModel)
* Slip model structure (e.g. SlipModel)
* Downsample parameter structure (e.g. SmoothParam)
  
A lot of the functions are designed as `Name-Value` pair inputs. Please refer to the following to see more.

---
 ### 1. okLoadData.m  
 * #### To read in GeoTiff files into Matlab environment using built-in function `readgeoraster()`

It is quite problematic to read ISCE product directly into Matlab. Therefore, convert ISCE product into Gdal GeoTiff with GDAL tools `gdal_translate` before reading in.  
  
Please let me know if there is an simpler way to do this. I appreciate any help!
```bash
gdal_translate -b 2 filt_topophase.unw.geo UnwrapPhase.tif
gdal_translate -b 1 los.rdr.geo Azimuth.tif
gdal_translate -b 2 los.rdr.geo Incidence.tif
gdal_translate -b 2 topophase.cor.geo Coherence.tif
```
After that, you can simply read them with the following. Note that **Coherence** is optional if you don't wish to mask your data with Coherence.  
#### Name-value pairs:  
* 'Displacement': ***Character or matrix.*** Input the displacement GeoTiff filename or an already read-in matrix
* 'Azimuth': ***Character or matrix.*** Input the azimuth angle GeoTiff filename or an already read-in matrix
* 'Incidence': ***Character or matrix.*** Input the incidence angle GeoTiff filename or an already read-in matrix
* 'Coherence': ***Character or matrix.*** Input the coherence GeoTiff filename or an already read-in matrix (This is optional)
* 'Processor': ***Character.*** Input the InSAR processor. (Right now only supports ISCE)
* 'LookDir': ***Character.*** Input the satellite/airplane looking direciton. (This is to convert ISCE azimuth angle convention to mine) ['r' or 'l']
```matlab
DataStruct = okLoadData('Displacement','UnwrapPhase.tif', ...
    'Azimuth','Azimuth.tif', ...
    'Incidence','Incidence.tif', ...
    'Coherence','Coherence.tif', ...
    'Processor','ISCE', ...
    'LookDir','r');
```
