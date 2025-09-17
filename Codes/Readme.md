#### Just want to clarify again here, `okada85.m` is a creation from [François Beauducel](https://github.com/IPGP/deformation-lib/tree/master/okada).  
#### My Matlab codes are designed to leverage its functionality.  

---
All Matlab variables are the format of ***structure***. It should be rather straight-forward to modify datasets by replacing data with the correct field.  
4 main data structures are used and created throughout the slip inversion:  
* Displacement structure (e.g. `DataStruct`)
* Fault model structure (e.g. `FaultModel`)
* Slip model structure (e.g. `SlipModel`)
* Downsample parameter structure (e.g. `SmoothParam`)
  
A lot of the functions are designed as `Name-Value` pair inputs. Please refer to the following to see more.  
Also, the structures are updated after each execution of the sub-routine, so set the ouput structure as the input as shown below.  

---
## 0. okPlot.m
It is highly recommended to use this plotting function to have a look of the stuff that is generated
#### Positional input:
* DataStruct: ***Structure.***: Input the sturcture that is going to be used to plot. 3 structures that are potentially used for this
  * Displacement structure (`DataStruct`)
  * Fault model structure (`FaultModel`)
  * Slip model structure (`SlipModel`)
* FigType: ***Character.*** Input what to plot. This is has to work with the input `DataStruct`. You cannot plot something that is not related to the data input  
  * '***Displ***': Plot displacement (`DataStruct`)
  * '***GNSS***': Plot GNSS displacement vectors (`DataStruct`)
  * '***Azimuth***': Plot azimuth anlge (`DataStruct`)
  * '***Incidence***': Plot incidence angle (`DataStruct`)
  * '***FaultModel***': Plot the fault model (`FaultModel`), use with sub-option `'displ'`
  * '***Dsample***': Plot the downsample result (`DataStruct`)
  * '***GreenFunc***': Plot the Okada Green's function (`FaultModel`), use with sub-options `'dset', 'displ'`
  * '***Smooth***': Plot the smoothing matrix (`FaultModel`)
  * '***Inversion***': Plot the inversion results (`SlipModel`), use with sub-options `'residual','slip1', 'slip2', 'totalslip', 'rake'`
  * '***AutoCorr***': Plot the auto-correlation result from InSAR (`DataStruct`), use with sub-options `'func','one'`
  * '***Covariance***': Plot the covariance matrix from InSAR auto-correlation (`DataStruct`)
#### Name-Value pairs:
* 'title': ***Character.***: Input figure title (default: ' ')
* 'clim': ***Vector.***: Input the colorbar limits
* 'cmap': ***Character.***: Input the colormap (default: 'default')
* 'displ': ***Structure or matrix***: Input the displacement, works with `FigType='FaultModel'`
* 'fnum': ***Numeric***: Input whether or not plot the patch number, works with `FigType='FaultModel'` (default: 0)
* 'markersize': ***Numeric.***: Input the size of the scatter plot (default: 5)
* 'dset': ***Character.***: Input the dataset to be plotted, works with `FigType='GreenFunc' & 'Inversion'` (default: 0)
* 'n': ***Numeric.***: Input which inversion result to plot, works with `FigType='Inversion'`
* 'residual': ***Numeric.***: Input whether or not plot the residual, works with `FigType='Inversion'` (default: 0)
* 'lcurve': ***Numeric.***: Input whether or not plot the L-curve, works with `FigType='Inversion'` (default: 0)
* 'slip1': ***Numeric.***: Input whether or not plot the first fault slip inversion, works with `FigType='Inversion'` (default: 0)
* 'slip2': ***Numeric.***: Input whether or not plot the second fault slip inversion, works with `FigType='Inversion'` (default: 0)
* 'totalslip': ***Numeric.***: Input whether or not plot the total fault slip inversion, works with `FigType='Inversion'` (default: 0)
* 'rake': ***Numeric.***: Input whether or not plot the rake inversion, works with `FigType='Inversion'` (default: 0)
* 'func': ***Character or cell of character.***: Input the fitted functions to be plotted, works with `FigType='AutoCorr'` (default: Plot all functions)
* 'one': ***Numeric.***: Input to plot the fitted functions in 1 figure or separate figures, works with `FigType='AutoCorr'` (default: 1)
* 'scale': ***Numeric***: Input the scaling factor for visualizing the GNSS vectors, works with `FigType='GNSS'` (default: 1)
* 'LineWidth': ***Numeric***: Input the line width of the GNSS vectors, works with `FigType='GNSS'` (default: 0.5)

### FigType = '***Displ***' || '***Azimuth***' || '***Incidence***'
```matlab
okPlot(DataStruct.Original,'Displ')
okPlot(DataStruct.Original,'Azimuth')
okPlot(DataStruct.Original,'Incidence')
```
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/f5da938c535b9843ede1a9afa722191e5b1374cb/Figure/okPlot_DisplAziInc.png)
### FigType = '***FaultModel***'
```matlab
okPlot(FaultModel,'FaultModel','displ',DataStruct.Subset)
```
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/f5da938c535b9843ede1a9afa722191e5b1374cb/Figure/okPlot_FaultModel.png)
### FigType = '***Dsample***'
```matlab
okPlot(DataStruct,'Dsample','MarkerSize',20)
```
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/f5da938c535b9843ede1a9afa722191e5b1374cb/Figure/okPlot_Dsample.png)
### FigType = '***GreenFunc***'
```matlab
okPlot(FaultModel,'GreenFunc','dset','GreenSS', ...
    'displ',DataStruct.Dsample, ...
    'MarkerSize',20)
```
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/f5da938c535b9843ede1a9afa722191e5b1374cb/Figure/okPlot_GreenFunc.png)
### FigType = '***Smooth***'
```matlab
okPlot(FaultModel,'Smooth')
```
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/5dd3cffa3a7a4b8734354a4c1db6ed5e5acdc1e6/Figure/okPlot_Smooth.png)
### FigType = '***Inversion***'
If `'n'` is not specified, the function will plot a several inversion results to give a quick skim
```matlab
okPlot(SlipModel,'Inversion','slip1',1,'n',5,'title','Dip-slip')
okPlot(SlipModel,'Inversion','slip2',1,'n',5,'title','Strike-slip')
okPlot(SlipModel,'Inversion','totalslip',1,'n',5)
okPlot(SlipModel,'Inversion','residual',1,'n',5,'MarkerSize',30)
okPlot(SlipModel,'Inversion','lcurve',1,'n',5,'title','L-curve')
```
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/5dd3cffa3a7a4b8734354a4c1db6ed5e5acdc1e6/Figure/okPlot_Inversion.png)
### FigType = '***AutoCorr***'
```matlab
okPlot(DataStruct.Subset,'AutoCorr','one',1)
```
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/9177e841dfaa50225c5279fe36eef7d416bc8b99/Figure/okPlot_AutoCorr.png)
### FigType = '***Covariance***'
```matlab
okPlot(DataStruct.Subset,'Covariance')
```
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/9177e841dfaa50225c5279fe36eef7d416bc8b99/Figure/okPlot_Covariance.png)

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
This is designed this way if users have their own masks wanted to apply  
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

### 5. okInSARCov.m
* ### To perform InSAR auto-correlation to produce covaraince matrix for weighting data
#### Positional input:
* DataStruct: ***Structure.*** Input DataStruct from `okLoadData.m`
* Dataset: ***Character.*** Input the field that contains the read-in displacement, azimuth ...
#### Name-Value pairs:
* 'deramp': ***Numeric.*** Input whether or not deramp first before auto-correlation (0 for no, default=1)
* 'func': ***Character or cell of characters.*** Input the functions used to fit the auto-correlation result (default='exp1','poly3','poly5','poly7','poly9')
```matlab
DataStruct = okInSARCov(DataStruct,'Subset', ...
    'deramp',1);
```

### 6. Down-sampling data 
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
### 7. okMakeInSARCovMat.m
* ### To make the covariance matrix from auto-correlation result for weighting InSAR data 
#### Positional input:
* DataStruct: ***Structure.*** Input DataStruct from `okLoadData.m`
* Dataset: ***Character.*** Input the field that contains the read-in displacement, azimuth ...
* func: ***Character.*** Input which fitted function result to be used to produce covariance matrix
```matlab
DataStruct = okMakeInSARCovMat(DataStruct,'Dsample','exp1');
```

### 8.1 okMakeFaultModel.m
* #### To make the fault geometry
#### Positional input:
* StartXYZ: ***Vector.***: Input the coordinate of the starting point of the fault patch
* Length: ***Numeric.***: Input the length of the fault
* Width: ***Numeric.***: Input the width (down-dip) of the fault
* Strike: ***Numeric.***: Input the strike of the fault (Right-hand rule)
* Dip: ***Numeric.***: Input the dip of the fault (Right-hand rule)
* PatchStrike: ***Numeric.***: Input the along-strike patch count
* PatchDip: ***Numeric.***: Input the along-dip patch count
```matlab
StartXYZ = [359099, 4253530, 0];
Length = 60000;
Width = 20000;
Strike = 242;
Dip = 80;
PatchStrike = 24;
PatchDip = 8;
FaultModel = okMakeFaultModel(StartXYZ,Length,Width,Strike,Dip,PatchStrike,PatchDip);
```

### 8.2 okCombineFaultModel.m
* #### To combine fault models with different patch sizes
#### Positional input:
* FaultModels: ***Cell.***: Input the fault models created from `okMakeFaultModel.m`
* Direction: ***String***: Input the direction of combining two faults
```matlab
Direction = 'strike';
Combined = okCombineFaultModel({FaultModel1,FaultModel2},Direction);

okPlot(Combined,'FaultModel','fnum',1);
```
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/c644b6744d2f572133ddb7de82b07eeb7bbfa4be/Figure/okCombineFaultModel.png)
---

### 9. okMakeGreenFunc.m
* #### To make the Okada Green's function  
Please be advised that the setting of the Green's function has a impact on fault slip inversion if using `Non-negative least squares` (as it forces all inversion result to be positive)  
#### Positional input:
* DataStruct: ***Structure.***: Input DataStruct
* Dataset: ***Character.***: Input the field that contains the read-in displacement, azimuth ...
* Direction: ***Character.***: Input which direction the Green's function is converted to
  * 'LOS': Convert to LOS direction based on the geometry in `DataStruct`
  * 'Azi': Convert to azimuth direction based on the geometry in `DataStruct`
  * '3D': Remain the Green's function in EW, NS and UD
* FaultModel: ***Structure.***: Input the FaultModel
* Rake: ***Numeric.***: Input the rake angle
  * Left-lateral: 0
  * Right-lateral: 180
  * Thrust: 90
  * Normal: -90
* Slip: ***Numeric.***: Input how much it slips (It can be zero to assume no slip in this rake direction)
* Opening: ***Numeric.***: Input how much tensile opening it is (It can be zero to assume no opening as a common practice in slip inversion)
* GreenFuncName. ***Character.***: Input the field name the Green's function is going to be stored in
```matlab
% Dip-slip
RakeDS = -90;
SlipDS = 1;
FaultModel = okMakeGreenFunc(DataStruct,'Dsample','LOS',FaultModel,RakeDS,SlipDS,0,'GreenDS');

% Strike-slip
RakeSS = 0;
SlipSS = 1;
FaultModel = okMakeGreenFunc(DataStruct,'Dsample','LOS',FaultModel,RakeSS,SlipSS,0,'GreenSS');
```
---

### 10. okMakeSmoothMat.m
* #### To make the smoothing matrix  
#### Positional input:
* FaultModel: ***Structure.***: Input FaultModel
#### Name-Value pairs:
* 'method': ***Character.***: Input the smoothing criterion
  * equidist: Treating all patches the same sizes and constant smoothing without considering the distance between them
  * dist-weighted: Smooth connected patches and weight them with their centroid distances
  * dist-based: Smooth all patches within the distance threshould and weight them with their centroid distances
* 'dist': ***Numeric.***: Input the distance threshold (Use with `dist-based` option)
```matlab
FaultModel = okMakeSmoothMat(FaultModel,'method','dist-weighted');
```
---

### 11. okInvertSlip.m
* #### To invert for the fault slip
To correctly use this function, you will need a very basic linear algebra.  
See below for different function input for different usages:
#### Positional input:
* DataStruct: ***Structure.***: Input DataStruct
* Dataset: ***Character.***: Input the field that contains the read-in displacement, azimuth ...
* FaultModel: ***Structure.***: Input the FaultModel
* GreenFuncDataset: ***Character.***: Input the corrsponding Green's function's field name
* GreenFuncPosition: ***Vector.***: Input where the Green's functions should be placed in the design matrix
* SmoothMatDataset: ***Character.*** Input the smoothing matrix field name stored in `FaultModel` (This is to give freedom for users using different smoother other than 2D Laplacian)
#### Name-Value pairs:
* 'solver': ***Character.***: Input the inversion solver
  * lsq: Ordinary Least Squares (default)
  * nnlsq: Non-negative Least Squares
  * fnnlsq: Fast non-negative Least Squares (From Bill Whiten (2025). nnls - Non negative least squares (https://www.mathworks.com/matlabcentral/fileexchange/38003-nnls-non-negative-least-squares), MATLAB Central File Exchange. Retrieved September 17, 2025.)
* 'smoothsearch': ***Vector.***: Input the search range of the smoothing constant (Leave blank to use the default one)
* 'weight': ***Character or cell of character or matrix.***: Input the weighting matrix. It can be the covariance matrix made in `okMakeInSARCovMat.m` or user defined weighting matrix (Leave blank to assume uniform weighting)
* 'lb': ***Vector.***: Input the lower bounds of the solutions
* 'ub': ***Vector.***: Input the upper bounds of the solutions
  
### Case 1: 1 LOS displacement, 2 Green's functions, InSAR covaraince matrix as data weights
This is a very common case for InSAR users to invert for fault slip. Here's how the matrices should be placed and the way this function is designed:  
  
<p align="center">
$$LOS$$ is the displacement </p>
<p align="center">
$$S$$ is the smoothing matrix </p>
<p align="center">
$$C$$ is the covariance matrix </p>
<p align="center">
$$G_{DS}$$ is the dip-slip Green's function </p>
<p align="center">
$$G_{SS}$$ is the strike-slip Green's function </p>  

And we are solving a linear system such that:
<p align="center">
$$d = Am$$ </p>  

The above linear system is weighted by $C^{-1}$ and smoothed by $S$:

$$(A^TC^{-1}A + \lambda^2S^TS)m = A^TC^{-1}d$$

The linear system can be written as an augmented format (convenient for coding):  
<p align="center">
$$\begin{bmatrix} C^{-1}G_{DS} & C^{-1}G_{SS} & 1 \\\ \lambda S & \lambda S & 0 \end{bmatrix}m = \begin{bmatrix} C^{-1}LOS \\\ 0 \end{bmatrix}$$

Here, the Green's functions are placed at matrix index `[11,12]`:  
<p align="center">
$$A_{index} = \begin{bmatrix} 11 & 12 & 13 \\\ 21 & 22 & 23 \end{bmatrix}$$  </p>

`S` will be automatically extended in `okInvertSlip,m` so don't worry about it. Therefore, your funciton input should be:  

```matlab
SmoothParam = [0.00001,0.0001,0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1,10];
GreenFuncPosition = [11,12];
SlipModel1 = okInvertSlip(DataStruct,'Dsample',FaultModel,{'GreenDS','GreenSS'},GreenFuncPosition,'SmoothMat', ...
    'solver','nnlsq', ...
    'smoothsearch',SmoothParam, ...
    'weight','Covariance');
```

### Case 2: 2 LOS displacements, 2 Green's functions
This is also a very common case for InSAR users to invert for fault slip where we have both ascending and descending observations
It is very similiar as the above
<p align="center">
$$LOS_{Asc}$$ is the ascending displacement </p>
<p align="center">
$$LOS_{Des}$$ is the descending displacement </p>
<p align="center">
$$S$$ is the smoothing matrix </p>
<p align="center">
$$C$$ is the covariance matrix </p>
<p align="center">
$$G_{AscDS}$$ is the ascending dip-slip Green's function </p>
<p align="center">
$$G_{AscSS}$$ is the ascending strike-slip Green's function </p>  
<p align="center">
$$G_{DesDS}$$ is the descending dip-slip Green's function </p>
<p align="center">
$$G_{DesSS}$$ is the descending strike-slip Green's function </p>  

The linear system is placed as such:
<p align="center">
$$\begin{bmatrix} C^{-1}G_{AscDS} & C^{-1}G_{AscSS} & 1 \\\ C^{-1}G_{DesDS} & C^{-1}G_{DesSS} & 1 \\\ \lambda S & \lambda S & 0 \end{bmatrix}m = \begin{bmatrix} C^{-1}LOS_{Asc} \\\ C^{-1}LOS_{Des} \\\ 0 \end{bmatrix}$$

Here, the Green's functions are placed at matrix index `[11,12,21,22]` with corrspondence with the displacement `[11,21]`: 
<p align="center">
$$d_{index} = \begin{bmatrix} 11 \\\ 21 \\\ 31 \end{bmatrix}$$ $$A_{index} = \begin{bmatrix} 11 & 12 & 13 \\\ 21 & 22 & 23 \\\ 31 & 32 & 33 \end{bmatrix}$$  </p>
Therefore, the function input will be:  

```matlab
SmoothParam = [0.00001,0.0001,0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1,10];
GreenFuncPosition = [11,12,21,22];
SlipModel = okInvertSlip({DataStructAsc,DataStructDes},{'Dsample','Dsample'}, ...
    FaultModel,{'GreenAscDS','GreenAscSS','GreenDesDS','GreenDesSS'},[11,12,21,22],'SmoothMat', ...
    'solver','nnlsq', ...
    'smoothsearch',SmoothParam, ...
    'weight',{'Covariance','Covariance'});
```

---

### 11. okOutputGMT.m
* #### To output result in GMT plottable format along with a slip inversion report
#### Positional input:
* SlipModel: ***Structure.***: Input SlipModel
* OutType: ***Character.***: Input the thing you want to output
  * 'slip': The slip inversion. 4 files asscociated with this will be made:
    * XXXX_Solution.txt: The inversion report
    * XXXX_Slip1_GMT.txt: GMT plottable fault slip associated with the first Green's function input
    * XXXX_Slip2_GMT.txt: GMT plottable fault slip associated with the second Green's function input
    * XXXX_TotalSlip_GMT.txt: GMT plottable total fault slip
  * 'lcurve': The L-curve search result
  * 'displ': The displacement observation and forward results
  * 'green': The Okada Green's function
#### Name-Value pairs:
* 'n': ***Numeric.***: Choose the slip inversion result (Based on the L-curve search)
* 'o': ***Character.*** The ouput file name prefix. (Leave blank the prefix will be `okSlipResult`)

```matlab
okOutputGMT(SlipModel,'slip', ...
    'n',5);
okOutputGMT(SlipModel,'lcurve', ...
    'n',5);
okOutputGMT(SlipModel,'green', ...
    'n',5);
okOutputGMT(SlipModel,'displ', ...
    'n',5);
```




