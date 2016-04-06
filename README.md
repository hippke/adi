# Accretion Disk Inversion
We are currently preparing a paper which, amongst other things, discusses the objects that cause dips in the light curve of a star with an accretion disk. This repository provides code and instructions to "inverse" the light curve, i.e. estimate the shape of the occulters. The code can read Kepler (e.g. K2-SFF) time-flux data and output light curves and videos.

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

### Inversion
Given the clean data, we can now perform the inversion. For each data point, we calculate the area required for a trial impact parameter. The limb darkening law is arbitrary; here we use the quadratic law:
```python
def quadratic_limb_darkening(impact, l1, l2):
    """Quadratic limb darkening. Kopal 1950, Harvard Col. Obs. Circ., 454, 1"""
    return 1 - l1 * (1 - impact) - l2 * (1 - impact) ** 2
```
The whole procedure will be done numerically, with a star including limb darkening being projected on a 2D numpy grid:
```python
def numerical_star(pixels, l1, l2):
    """Provides a numerical array of flux-values for a limb-darkened star"""
    grid = numpy.zeros([2 * pixels, 2 * pixels])
    for width in range(2 * pixels):
        for height in range(2 * pixels):
            distance_total = math.sqrt(((pixels - width) ** 2) +
                ((pixels - height) ** 2))
            if distance_total < pixels:
                circle_position = distance_total / float(pixels)
                grid[width, height] = quadratic_limb_darkening(
                    1 - circle_position, l1, l2)
    return grid
```

We inject the trial occultar array into the star grid:
```python
def star_with_occult(pixels, l1, l2, occultarray):
    """Fetches limb-darkened star and overplots occulter percentages"""
    star_grid = numpy.copy(buffered_grid)  # numerical_star(pixels, l1, l2)
    for x_value in range(len(occultarray)):
        height = int(occultarray[x_value])
        # Erase area under curve
        for y_value in range(height, pixels * 2):
            star_grid[y_value, x_value] = 0.  # Set opacity [0..1] here
    return star_grid
```

And can obtain an image (here with a colormap for illustration purpose only):


![Image](http://www.jaekle.info/c1.png "Img1")

![Image](http://www.jaekle.info/c3.png "Img1")
