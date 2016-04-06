# Accretion Disk Inversion
We are currently preparing a paper which, amongst other things, discusses the objects that cause dips in the light curve of a star with an accretion disk. This repository provides code and instructions to "inverse" the light curve, i.e. estimate the shape of the occulters. The code can read Kepler (e.g. K2-SFF) time-flux data and output lightcurvers and videos.

## Tutorial
### Read K2SFF data and generate spline fit
Some smoothing is useful as we need a continous lightcurve, and don't want to fit out noise features or spots. Here we use a spline, but of course any other method can be used, e.g. a sliding median.
```python
time, flux = read_k2sff_file(filename)
spline_time, spline_flux = make_spline_fit(
    time, flux, allowed_errors, upsample_factor)
```
We can plot the resulting scatter data (red symbols) and overlay the spline fit:
![Image](http://www.jaekle.info/c2.png "Img1")


![Image](http://www.jaekle.info/c1.png "Img1")

![Image](http://www.jaekle.info/c3.png "Img1")
