The data could be download from my [Google folder](https://drive.google.com/drive/folders/1XrQMD9pLxDESoMZopdGDDkvZHXZlDM0E?usp=sharing)  

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
