import numpy
import math
import matplotlib
import scipy.interpolate as inter
import matplotlib.pyplot as plt
import matplotlib.colors as col


class star:
    """Sets and gets ld'd scale height for rectangle corresponding to flux"""
    def __init__(self, ld_steps, l1, l2):
        buffered_rectangle_inverse = numpy.array([])
        lookup_grid_ld = numerical_star(ld_steps, l1, l2)
        lookup_allflux = numpy.sum(lookup_grid_ld)
        for i in range(ld_steps + 1):
            impact = i / ld_steps
            remaining_flux = \
                1 - (numpy.sum(lookup_grid_ld[:,:i * 2]) / lookup_allflux)
            buffered_rectangle_inverse = numpy.append(
                buffered_rectangle_inverse, remaining_flux)
        self.__x = buffered_rectangle_inverse

    def inverse(self, ld_steps, required_total_flux):
        """Returns scale height of rectangle corresponding to req flux value"""
        current_column_occult_ratio = 0
        for i in range(ld_steps + 1):
            if required_total_flux > self.__x[i]:
                current_column_occult_ratio = self.__x[i]
                break
        return current_column_occult_ratio


def quadratic_limb_darkening(impact, l1, l2):
    """Quadratic limb darkening. Kopal 1950, Harvard Col. Obs. Circ., 454, 1"""
    impact = math.cos(impact)
    return 1 - l1 * (1 - impact) - l2 * (1 - impact) ** 2


def numerical_star(pixels, l1, l2):
    """Provides a numerical array of flux-values for a limb-darkened star"""
    grid = numpy.zeros([2 * pixels, 2 * pixels])
    for width in range(2 * pixels):
        for height in range(2 * pixels):
            distance_total = \
                math.sqrt(((pixels - width) ** 2) + ((pixels - height) ** 2))
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
    solution = 1.
    for i in range(steps + 1):
        current_ld = quadratic_limb_darkening(i / steps, l1, l2)
        # Speed increase if list is created once and then buffered
        # print(i / steps, current_ld)
        if current_ld < value_to_find:
            solution = i / steps
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
#filename = 'hlsp_k2sff_k2_lightcurve_204757338-c02_kepler_v1_llc-default-aper.txt'  # first test
#filename = 'hlsp_k2sff_k2_lightcurve_203343161-c02_kepler_v1_llc-default-aper.txt'
filename = 'hlsp_k2sff_k2_lightcurve_204137184-c02_kepler_v1_llc-default-aper.txt'
upsample_factor = 10
allowed_errors = 0.005  # Check and adjust spline fit visually
allowed_errors = 0.01
stellar_radius = 25  # Unit: pixels
ld_steps = 250  # Numerical steps for limb darkening inversion
l1 = 0.4364  # Quadratic limb darkening, parameter 1
l2 = 0.3467  # Quadratic limb darkening, parameter 2
alpha = 2.

# Speed tricks
buffered_grid = numerical_star(stellar_radius, l1, l2)
print('Generating ld star...')
buffered_star = star(ld_steps, l1, l2)
lookup_grid_ld = numerical_star(ld_steps, l1, l2)
lookup_allflux = numpy.sum(lookup_grid_ld)
required_total_flux = 1.  # Just an initial guess
occulted_flux_ratio = 0.  # Just an initial guess
unocculted_flux = numpy.sum(buffered_grid)  # Total flux without occultation

# Colormap with black background and red-yellow star
startcolor = 'black'
midcolor = 'yellow'
endcolor = (1, 0.5, 0)  # yellow
cmap = col.LinearSegmentedColormap.from_list(
    'own2', [startcolor,midcolor, midcolor, midcolor, endcolor])

# Read K2SFF data and generate spline fit
time, flux = read_k2sff_file(filename)
# time = time[2200:]
# flux = flux[2200:]
print('Performing spline fit...')
spline_time, spline_flux = make_spline_fit(
    time, flux, allowed_errors, upsample_factor)

# Check fit visually
plt.scatter(time, flux, color='r', alpha=0.5, s=20)
plt.plot(spline_time, spline_flux, 'k', linewidth=2.)
plt.plot((0, 5000), (1, 1), 'k', linewidth=0.5)
plt.xlim(numpy.amin(time), numpy.amax(time))
# plt.ylim(numpy.amin(flux), numpy.amax(flux))
plt.ylim(numpy.amin(flux), 1.1)
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
# plt.show()
plt.clf()

# Create inversion
occultarray = numpy.empty(stellar_radius * 2)  # Start with empty array
occultarray.fill(stellar_radius * 2)  # Populate with ld star background
current_column_occult_ratio = 0.
list_wanted = []
list_fitted = []

print('Performing inversion...')

for currentvalue in range(len(spline_flux)):
    # Look ahead: in one radius distance, what occultation level will occur?
    if currentvalue < stellar_radius:
        required_total_flux = spline_flux[currentvalue]
    else:
        required_total_flux = spline_flux[currentvalue - stellar_radius]

    current_column_occult_ratio = buffered_star.inverse(ld_steps, required_total_flux)   
    current_column_occult_ratio = current_column_occult_ratio * alpha + (1 - alpha)
    # Assure value is within 0..1
    current_column_occult_ratio = max(0, min(current_column_occult_ratio, 1))

    # Insert best bet ratio in first position of new array
    occultarray[0] = (current_column_occult_ratio) * stellar_radius * 2

    # Produce a grid with the new occultation sitation
    occultedflux = star_with_occult(stellar_radius, l1, l2, occultarray)
    occulted_flux_ratio = numpy.sum(occultedflux) / unocculted_flux

    # Roll array to the right
    occultarray = numpy.roll(occultarray, 1)

    # If you want to create a video, then switch on saving the frames
    # if currentvalue % 20 == 0:
        # save_videoframe(occultedflux, currentvalue)  # Slow!

    if currentvalue % 1000 == 0:  # Speedup when not printing everything
        print('progress', "{0:.1f}%".format(currentvalue / len(spline_flux) * 100))
        # print(currentvalue, required_total_flux, occulted_flux_ratio,
            # current_column_occult_ratio)
    list_wanted.append(str(required_total_flux))
    list_fitted.append(str(occulted_flux_ratio))

print('Inversion complete. Making figure...')

# Make plot to compare data and fit
plt.scatter(time, flux, color='r', alpha=0.5, s=20)
print_offset = (2 * stellar_radius / len(spline_time)) * \
    (numpy.amax(spline_time) - numpy.amin(spline_time))
plt.plot(spline_time - print_offset, list_fitted, 'k', linewidth=2.)
plt.plot((0, 5000), (1, 1), 'k', linewidth=0.5)
plt.xlim(numpy.amin(time), numpy.amax(time))
plt.ylim(0.5, 1.1)
# plt.ylim(0.6, 1.05)
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

"""
# Make zoom plot into last part of dataset
plt.scatter(time, flux, color='r', alpha=0.5, s=20)
print_offset = (2 * stellar_radius / len(spline_time)) * \
    (numpy.amax(spline_time) - numpy.amin(spline_time))
plt.plot(spline_time - print_offset, list_fitted, 'k', linewidth=2.)
plt.plot((0, 5000), (1, 1), 'k', linewidth=0.5)
plt.xlim(2124.5, numpy.amax(time))
plt.ylim(0.5, 1.1)
#plt.ylim(0.8, 1.05)
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
"""
