#### Just want to clarify again here, `okada85.m` is a creation from [okada85.m](https://github.com/IPGP/deformation-lib/tree/master/okada). My Matlab codes are designed to leverage its functionality.  

---
All Matlab variables are the format of ***structure***. It should be rather straight-forward to modify datasets by replacing data with the correct field.  
4 main data structures are used and created throughout the slip inversion:  
* Displacement structure (e.g. DataStruct)
* Fault model structure (e.g. FaultModel)
* Slip model structure (e.g. SlipModel)
* Downsample parameter structure (e.g. SmoothParam)
  
A lot of the functions are designed as `Name-Value` pair inputs. Please refer to the following to see more.  
Also, the structures are updated after each execution of the sub-routine, so set the ouput structure as the input as shown below.  

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
#### Name-Value pairs:  
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
---

### 2. okMaskData.m
* #### To mask data with coherence or a read-in matrix with the same size  
#### Positional input:  
* DataStruct: ***Structure.*** Input DataStruct from `okLoadData.m`
* Dataset: ***Character.*** Input the field that contains the read-in displacement, azimuth ...
* TobeMasked: ***Character.*** Input the field to be masked
* MaskCriterion: ***Matrix.*** Input the matrix to be used as the mask  
#### Name-Value pairs:
* 'Threshold': ***Numeric.*** Input the set threshold for `MaskCriterion` (Leave blank to assume `MaskCriterion` is already a mask)  
```matlab
Threshold = 0.35;
DataStruct = okMaskData(DataStruct,'Original','Displ',DataStruct.Original.Coherence, ...
    'Threshold',Threshold);
% Or if MaskCriterion is already a mask then:
DataStruct = okMaskData(DataStruct,'Original','Displ',Mask);
```
---

### 3. okMakeDataSubset.m
* #### To make a subset of the data strucuture (Cropping)  
#### Positional input:  
* DataStruct: ***Structure.*** Input DataStruct from `okLoadData.m`
* Dataset: ***Character.*** Input the field that contains the read-in displacement, azimuth ...
* rmin: ***Numeric.*** Input the minimum row index
* rmax: ***Numeric.*** Input the maximum row index
* cmin: ***Numeric.*** Input the minimum column index
* cmax: ***Numeric.*** Input the maximum column index
```matlab
rmin = 3500;
rmax = 5932;
cmin = 1;
cmax = 4926;
DataStruct = okMakeDataSubset(DataStruct,'Original',rmin,rmax,cmin,cmax);
```
---

### 4. okLLtoLocalXY.m
* #### To convert longitude and latitude to local XY coordinates given a longitude origin
#### Positional input:
* DataStruct: ***Structure.*** Input DataStruct
* Dataset: ***Character.*** Input the field that contains the read-in displacement, azimuth ...
* lon_origin: ***Numeric.*** Input the longitude origin for conversion
```matlab
% After okMakeDataSubset.m, there is a new field called 'Subset' that contains the subset of your data
lon_origin = 38;
DataStruct = okLLtoLocalXY(DataStruct,'Subset',lon_origin);
```
---

### 5. Down-sampling data 
To down-sample the data using quadtree algorithm.  
This will be a trial-and-error process until you find the best downsampled results for your case.  

- ### okMakeDsampleParam.m  
Make the parameter structure for `okDsample.m`

#### Positional input:
* FuncType: ***Character.*** Input the downsampling algorithm. (Right now only supports quadtree)
#### Name-Value pairs:
* 'Criterion': ***Character.*** Input the quadtree criterion ['variance' or 'curvature']
  * variance: Method used in Jonsson (2002)
  * curvature: Method used in Simons et al. (2002)
* 'Tolerance': ***Numeric.*** Input the variance threshold in each leaf
* 'LeafMinSize': ***Numeric.*** Input the minimum size 1 leaf can be
* 'LeafMaxSize': ***Numeric.*** Input the maximum size 1 leaf can be
* 'NaNPixelAllow': ***Numeric.*** Input the maximum NaN pixel ration 1 leaf can have

- ### okDsample.m
Calls the downsampling algorithm based on the parameters made in `okMakeDsampleParam.m`

#### Positional input:
* DataStruct: ***Structure.*** Input DataStruct
* Dataset: ***Character.*** Input the field that contains the read-in displacement, azimuth ...
* FuncParam: ***Structure.*** Input the parameter structure made in `okMakeDsampleParam.m`

- ### okQuadtreeD.m
Recursive function that operates the quadtree downsampling

```matlab
% Make the downsample parameter structure
DsampleParam = okMakeDsampleParam('quadtree', ...
    'Criterion','curvature', ...
    'Tolerance',0.00007, ...
    'LeafMinSize',1000, ...
    'LeafMaxSize',200000, ...
    'NaNPixelAllow',0.5);
% Downsample the data based on the selected algorithm 
DataStruct = okDsample(DataStruct,'Subset',DsampleParam);
```
---




