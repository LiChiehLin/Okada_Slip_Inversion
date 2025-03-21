# Okada_Slip_Inversion
This is a suite of Matlab routines for using InSAR data to invert for fault slip using Okada solution [(Okada, 1985)](https://pubs.geoscienceworld.org/ssa/bssa/article/75/4/1135/118782/Surface-deformation-due-to-shear-and-tensile).  
The fundamental Okada solution is from **François Beauducel** [okada85.m](https://github.com/IPGP/deformation-lib/tree/master/okada).  I DON'T own any credit of the creation of okada85.m.  

---
The basic workflow is as follows:  
  
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/7feebc821cd85102756997bff945e8ab1af74999/Figure/Workflow.png)

## Below is an example for [2020 Mw 6.8 Elazig earthquake](https://github.com/LiChiehLin/Okada_Slip_Inversion/tree/24275589c24899d169ea9e3841b004c4a88fe686/Example_elazig_SinglePlane)
1. Convert ISCE product to Matlab readable GeoTiff using GDAL tools:
```sh
gdal_translate -b 2 filt_topophase.unw.geo UnwrapPhase.tif
gdal_translate -b 1 los.rdr.geo Azimuth.tif
gdal_translate -b 2 los.rdr.geo Incidence.tif
gdal_translate -b 2 topophase.cor.geo Coherence.tif
```

2. Use `okLoadData.m` to read these files into Matlab environment
```matlab
DataStruct = okLoadData('Displacement','UnwrapPhase.tif', ...
    'Azimuth','Azimuth.tif', ...
    'Incidence','Incidence.tif', ...
    'Coherence','Coherence.tif', ...
    'Processor','ISCE', ...
    'LookDir','r');
```
Note that, since ISCE outputs radians rather than meters, you have to manually convert radian to meters by:  
<p align="center">
$$meters = \frac{φ*λ}{4π}$$  </p>
<p align="center">
where φ is displacement radians, λ is wavelength in meters </p>  

```matlab
wavel = 0.055;
DataStruct.Original.Displ = DataStruct.Original.Displ.*wavel./4./pi;
```

3. Use `okMaskData.m` to mask noisy pixels using coherence  
```matlab
Threshold = 0.35;
DataStruct = okMaskData(DataStruct,'Original','Displ',DataStruct.Original.Coherence,'Threshold',Threshold);
```

4. Use `okMakeDataSubset.m` to crop the data into a smaller piece  
```matlab
rmin = 3500;
rmax = 5932;
cmin = 1;
cmax = 4926;
DataStruct = okMakeDataSubset(DataStruct,'Original',rmin,rmax,cmin,cmax);
```

5. Use `okLLtoLocalXY.m` to convert Lon Lat to Local XY coordinate system
```matlab
lon_origin = 38;
DataStruct = okLLtoLocalXY(DataStruct,'Subset',lon_origin);
```

6. Use `okMakeDsampleParam.m` to make downsample parameters and use `okDsample.m` to downsample the InSAR displacement  
Currently there is only **Quadtree** downsampling algorithm and 2 downsampling criterion can be chosen from:  
* **'variance'**: [Jonsson (2002). Modeling volcano and earthquake deformation from satellite radar interferometric observations. Stanford University.](https://www.proquest.com/docview/305523554?pq-origsite=gscholar&fromopenview=true&sourcetype=Dissertations%20&%20Theses)    
* **'curvature'**: [Simons et al. (2002). Coseismic deformation from the 1999 M w 7.1 Hector Mine, California, earthquake as inferred from InSAR and GPS observations. Bulletin of the Seismological Society of America, 92(4), 1390-1402.](https://pubs.geoscienceworld.org/ssa/bssa/article/92/4/1390/120788/Coseismic-Deformation-from-the-1999-Mw-7-1-Hector)  
```matlab
DsampleParam = okMakeDsampleParam('quadtree', ...
    'Criterion','curvature', ...
    'Tolerance',0.00007, ...
    'LeafMinSize',1000, ...
    'LeafMaxSize',200000, ...
    'NaNPixelAllow',0.5);
DataStruct = okDsample(DataStruct,'Subset',DsampleParam);
```

7. Use `okMakefaultModel.m` to make the fault geometry
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

8. Use `okMakeGreenFunc.m` to make Strike-slip and Dip-slip Green's functions
```matlab
RakeDS = -90;
RakeSS = 0;
SlipDS = 1;
SlipSS = 1;
FaultModel = okMakeGreenFunc(DataStruct,'Dsample','LOS',FaultModel,RakeDS,SlipDS,0,'GreenDS');
FaultModel = okMakeGreenFunc(DataStruct,'Dsample','LOS',FaultModel,RakeSS,SlipSS,0,'GreenSS');
```

9. Use `okMakeSmoothMat.m` to make the smoothing matrix
Currently there is only the 2D Laplacian smoother  
```matlab
FaultModel = okMakeSmoothMat(FaultModel);
```

10. Use `okInvertSlip.m` to invert for fault slip
Two solvers to choose from:
* 'lsq': Ordinary Least Squares
* 'nnlsq': Non-negative Least Squares
```matlab
SmoothParam = [0.00001,0.0001,0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1,10];
SlipModel = okInvertSlip(DataStruct,'Dsample',FaultModel,{'GreenDS','GreenSS'},[11,12],'SmoothMat', ...
    'solver','nnlsq', ...
    'smoothsearch',SmoothParam);
```

11. Use `okOutputGMT.m` to output to GMT plottable format (psxy -L)
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

---

### Results:
#### L-curve:
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/68a1fd2796768bbd04747121af4c791c03af7305/Figure/Lcurve.png)

#### Fault slip
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/3961e537622757f0c652d0debde8539f9c841fda/Figure/TotalSlip.png)

#### Residual
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/fb2534304298bb9a75073b6fe4a66491950c619b/Figure/Residual.png)
