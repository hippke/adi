import numpy
import math
import csv
import matplotlib
import scipy.interpolate as inter
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def quadratic_limb_darkening(impact, l1, l2):
    """Quadratic limb darkening. Kopal 1950, Harvard Col. Obs. Circ., 454, 1"""
    impact = math.cos(1 - impact)
    return 1 - l1 * (1 - impact) - l2 * (1 - impact) ** 2


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


def star_with_occult(pixels, l1, l2, occultarray):
    """Fetches limb-darkened star and overplots occulter percentages"""
    star_grid = numpy.copy(buffered_grid)  # numerical_star(pixels, l1, l2)
    for x_value in range(len(occultarray)):
        height = int(occultarray[x_value])
        # Erase area under curve
        for y_value in range(height, pixels * 2):
            star_grid[y_value, x_value] = 0.  # Set opacity [0..1] here
    return star_grid


def get_ld_inversion(value_to_find, l1, l2, steps):
    """Returns the impact parameter for a given flux level"""
    for i in range(steps + 1):
        current_ld = quadratic_limb_darkening(1 - (i / steps), l1, l2)
        # Speed increase if list is created once and then buffered
        if current_ld < value_to_find:
            solution = 1 - (i / steps)
            break
    return solution


def read_k2sff_file(filename):
    """Reads K2SFF file from given name and returns time, flux numpy arrays"""
    time = numpy.loadtxt(filename, delimiter=',', skiprows=1, usecols=[0])
    flux = numpy.loadtxt(filename, delimiter=',', skiprows=1, usecols=[1])
    return time, flux


def make_spline_fit(time, flux, allowed_errors, upsample_factor):
    """Returns the spline to a given time, flux dataset"""
    spline = inter.UnivariateSpline(time, flux, s=allowed_errors)
    output_time = numpy.linspace(
        min(time), max(time), len(time) * upsample_factor)
    return output_time, spline(output_time)


def save_videoframe(fluxarray, framenumber):
    fig = plt.figure(frameon=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    fig.add_axes(ax)
    ax.imshow(fluxarray, aspect='auto', cmap=cmap)
    # for small sizes add ,interpolation='none'
    fig.savefig('fig' + str(framenumber) + '.png', dpi=100)
    plt.close()


"""Start of main routine"""
# Set parameters
filename = 'hlsp_k2sff_k2_lightcurve_204757338-c02_kepler_v1_llc-default-aper.txt'
upsample_factor = 50
allowed_errors = 0.005  # Check and adjust spline fit visually
stellar_radius = 100  # Unit: pixels
ld_steps = 1000  # Numerical steps for limb darkening inversion
l1 = 0.4364  # Quadratic limb darkening, parameter 1
l2 = 0.3467  # Quadratic limb darkening, parameter 2
# Curve morphology: We do not transit over the center of the star
# Include this value as a free parameter in the fit. Typical values are 0.95--1
# Can fit out the delta within the flux from deviating from a rectangle
correction_factor = 0.985
buffered_grid = numerical_star(stellar_radius, l1, l2)  # Optional speed trick
# Colormap with black background and red-yellow star
cmap = LinearSegmentedColormap.from_list('mycmap', [(0 / 3.0, 'black'),
    (1 / 3.0, 'red'), (3 / 3.0, 'yellow')])

# Read K2SFF data and generate spline fit
time, flux = read_k2sff_file(filename)
spline_time, spline_flux = make_spline_fit(
    time, flux, allowed_errors, upsample_factor)

# Check fit visually
plt.scatter(time, flux, color='r', alpha=0.5, s=20)
plt.plot(spline_time, spline_flux, 'k', linewidth=2.)
plt.plot((0, 5000), (1, 1), 'k', linewidth=0.5)
plt.xlim(numpy.amin(time), numpy.amax(time))
# plt.ylim(numpy.amin(flux), numpy.amax(flux))
plt.ylim(0.9, 1.05)
ax = plt.axes()
fmt = matplotlib.ticker.ScalarFormatter(useOffset=False)
fmt.set_scientific(False)
ax.xaxis.set_major_formatter(fmt)
ax.yaxis.set_major_formatter(fmt)
plt.tick_params(axis='both', which='major', labelsize=16)
ax.tick_params(direction='out')
plt.xlabel('Time (days)', fontsize=16, fontweight='bold')
plt.ylabel('Relative Flux', fontsize=16, fontweight='bold')
plt.savefig("fig_spline_scatter_" + filename + ".pdf", bbox_inches='tight')
plt.show()
plt.clf()
# If the fit is bad, adjust "allowed_errors"

# Create inversion
required_total_flux = 1.  # Just an initial guess
occulted_flux_ratio = 0.  # Just an initial guess
unocculted_flux = numpy.sum(buffered_grid)  # Total flux without occultation
occultarray = numpy.empty(stellar_radius * 2)  # Start with empty array
occultarray.fill(stellar_radius * 2)  # Populate with ld star background
current_column_occult_ratio = 0.
list_wanted = []
list_fitted = []

for currentvalue in range(len(spline_flux)):
    # Look ahead: in one radius distance, what occultation level will occur?
    if currentvalue < stellar_radius:
        required_total_flux = spline_flux[currentvalue]
    else:
        required_total_flux = spline_flux[currentvalue - stellar_radius]

    current_column_occult_ratio = 1 - get_ld_inversion(
        required_total_flux * correction_factor, l1, l2, ld_steps)

    # Make sure the occultation amount is in valid range
    if current_column_occult_ratio < 0:
        current_column_occult_ratio = 0
    if current_column_occult_ratio > 1:
        current_column_occult_ratio = 1.

    # Insert best bet ratio in first position of new array
    occultarray[0] = (1 - current_column_occult_ratio) * stellar_radius * 2

    # Produce a grid with the new occultation sitation
    occultedflux = star_with_occult(stellar_radius, l1, l2, occultarray)
    occulted_flux_ratio = numpy.sum(occultedflux) / unocculted_flux

    # Roll array to the right
    occultarray = numpy.roll(occultarray, 1)

    # If you want to create a video, then switch on saving the frames
    #if spline_time[currentvalue] > 2124.5 and currentvalue % 20 == 0:
        #save_videoframe(occultedflux, currentvalue)  # Slow!

    if currentvalue % 25 == 0:  # Speedup when not printing everything
        print(currentvalue, required_total_flux, occulted_flux_ratio,
            current_column_occult_ratio)
    list_wanted.append(str(required_total_flux))
    list_fitted.append(str(occulted_flux_ratio))

# Make plot to compare data and fit
plt.scatter(time, flux, color='r', alpha=0.5, s=20)
print_offset = (2 * stellar_radius / len(spline_time)) * \
    (numpy.amax(spline_time) - numpy.amin(spline_time))
plt.plot(spline_time - print_offset, list_fitted, 'k', linewidth=2.)
plt.plot((0, 5000), (1, 1), 'k', linewidth=0.5)
plt.xlim(numpy.amin(time), numpy.amax(time))
# plt.ylim(numpy.amin(flux), numpy.amax(flux))
plt.ylim(0.9, 1.05)
ax = plt.axes()
fmt = matplotlib.ticker.ScalarFormatter(useOffset=False)
fmt.set_scientific(False)
ax.xaxis.set_major_formatter(fmt)
ax.yaxis.set_major_formatter(fmt)
plt.tick_params(axis='both', which='major', labelsize=16)
ax.tick_params(direction='out')
plt.xlabel('Time (days)', fontsize=16, fontweight='bold')
plt.ylabel('Relative Flux', fontsize=16, fontweight='bold')
plt.savefig("fig_inversion_" + filename + ".pdf", bbox_inches='tight')
plt.show()
plt.clf()

# Make zoom plot into last part of dataset
plt.scatter(time, flux, color='r', alpha=0.5, s=20)
print_offset = (2 * stellar_radius / len(spline_time)) * \
    (numpy.amax(spline_time) - numpy.amin(spline_time))
plt.plot(spline_time - print_offset, list_fitted, 'k', linewidth=2.)
plt.plot((0, 5000), (1, 1), 'k', linewidth=0.5)
plt.xlim(2124.5, numpy.amax(time))
# plt.ylim(numpy.amin(flux), numpy.amax(flux))
plt.ylim(0.9, 1.05)
ax = plt.axes()
fmt = matplotlib.ticker.ScalarFormatter(useOffset=False)
fmt.set_scientific(False)
ax.xaxis.set_major_formatter(fmt)
ax.yaxis.set_major_formatter(fmt)
plt.tick_params(axis='both', which='major', labelsize=16)
ax.tick_params(direction='out')
plt.xlabel('Time (days)', fontsize=16, fontweight='bold')
plt.ylabel('Relative Flux', fontsize=16, fontweight='bold')
plt.savefig("fig_inversion_zoom_" + filename + ".pdf", bbox_inches='tight')
plt.show()
plt.clf()
