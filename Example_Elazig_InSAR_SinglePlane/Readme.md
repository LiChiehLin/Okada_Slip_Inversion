The data could be downloaded from my [Google folder](https://drive.google.com/drive/folders/1XrQMD9pLxDESoMZopdGDDkvZHXZlDM0E?usp=sharing)   
Download the files and the [Codes/](https://github.com/LiChiehLin/Okada_Slip_Inversion/tree/main/Codes) and arrange them as this:

```bash
├── function/
│   ├── okada85.m
│   ├── okCombineFaultModel.m
│   ├── okDsample.m
│   ├── okInvertSlip.m
│   ├── okLLtoLocalXY.m
│   ├── okLoadData.m
│   ├── okLocalXYtoLL.m
│   ├── okMakeDataSubset.m
│   ├── okMakeDsampleParam.m
│   ├── okMakeFaultModel.m
│   ├── okMakeGreenFunc.m
│   ├── okMakeSmoothMat.m
│   ├── okMaskData.m
│   ├── okPlot
│   └── okQuadtreeD.m
├── Azimuth.tif
├── Coherence.tif
├── Incidence.tif
├── UnwrapPhase.tif
└── main_elazig_example.m
```

Add the routines to your working environment path:
```matlab
addpath(genpath('function/'))
```
You may run `main_elazig_example.m` seciton by section or all the way through.  
It should work all the way through and users can modify this main script to fit with their own needs.  

---
Use `okPlot.m` to assist you determine where to make the fault geometry. Do this after you convert to local coordinate system

```matlab
okPlot(DataStruct.Subset,'Displ')
```
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/665824cfa5cd65f14e2c4239acae3504818aaf7b/Figure/Displacement.png)

To visualize the quadtree downsampled result
```matlab
okPlot(DataStruct,'Dsample','MarkerSize',20)
```
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/665824cfa5cd65f14e2c4239acae3504818aaf7b/Figure/Quadtree.png)
