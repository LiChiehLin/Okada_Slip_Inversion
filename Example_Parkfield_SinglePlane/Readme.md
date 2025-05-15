## Example file for using GNSS data to invert for the fault slip distribution
- EQ: M6.0 2004 Parkfield earthquake  
- Data: ***[Houlié, Nicolas, Douglas Dreger, and Ahyi Kim. "GPS source solution of the 2004 Parkfield earthquake." Scientific reports 4.1 (2014): 3646.](https://www.nature.com/articles/srep03646)***

---
In the example file, I set all **vertical displacement** to `nan` and **did not** solve for the dip-slip component.  
You are encouraged to play with it and see the different results.

---
### Results:

#### L-curve
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/e3ae4a377dafdf582bc5a4b1f614f62f20e7ccc3/Figure/Parkfield_Lcurve.png)

#### Fault slip distribution
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/e3ae4a377dafdf582bc5a4b1f614f62f20e7ccc3/Figure/Parkfield_faultslip.png)

#### Model fit
![Example](https://github.com/LiChiehLin/Okada_Slip_Inversion/blob/e3ae4a377dafdf582bc5a4b1f614f62f20e7ccc3/Figure/Parkfield_modelfit.png)

---
### GMT plotting
You may use GMT to plot the fault slip distribution with `gmt psxy -L -C` command. Something like this:
```shell
gmt makecpt -Chot -T0/0.6/0.01 -I -Z > slip.cpt
gmt psbasemap -R-121/-120/35.5/36.2 -JM12c -K -B0.5 > test.ps
gmt psxy Parkfield_Slip1_GMT.txt -R -J -O -K -Cslip.cpt -L -t20 >> test.ps

echo | gmt psxy -R -J -O >> test.ps
gmt psconvert -Tg -A -Z test.ps
```
However, you won't be able to see anything because the dip angle was set to 90 degrees. You can change to different dip angles and see the effect.
