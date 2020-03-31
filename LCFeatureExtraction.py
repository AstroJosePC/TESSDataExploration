import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve, convolve_fft
from lightkurve import LightCurve

default_stddev, default_gauss_x = 3, 37


def smooth_lightcurve(orig_lc: LightCurve, kernel=None, gauss_size=None, gauss_std=None,
                      break_tolerance=5, convolution='convolve', **kwargs):
    """
    Smooth a light curve with a given kernel through convolution.

    :param orig_lc: light curve to be smoothed.
    :param kernel: The convolution kernel. Must be 1-Dimensional. Default Astropy's Gaussian1DKernel
    :param gauss_size: Size of the default Gaussian kernel. Default: SQRT[ light curve length ]
    :param gauss_std: Standard deviation of the default Gaussian kernel. Default: gauss_size / 10.0
    :param break_tolerance: tolerance used to break light curve into sub-lightcurves
    :param convolution: convolution function; default astropy.stats..convolution
    :param kwargs: Keyword arguments to pass to the convolution function
    :return:
    """
    if kernel is None:
        # Create Gaussian kernel
        if gauss_size is None:
            # Define odd gauss size ; method found from another project
            gauss_size = int(np.sqrt(orig_lc.flux.size))
            if gauss_size % 2 == 0:
                gauss_size += 1
        if gauss_std is None:
            # Define standard deviation of Gaussian model; experimentally determined
            gauss_std = gauss_size / 10.0
        kernel = Gaussian1DKernel(stddev=gauss_std, x_size=gauss_size)

    # Choose convolution method
    if convolution == 'convolve':
        convolve_teq = convolve
    elif convolution == 'fft':
        convolve_teq = convolve_fft
    else:
        raise ValueError(f'Wrong convolution technique "{convolution}"')

    # Split the lightcurve into segments by finding large gaps in time
    # This segmentation method is borrowed from lightkurve flatten method
    dt = orig_lc.time[1:] - orig_lc.time[0:-1]
    cut = np.where(dt > break_tolerance * np.nanmedian(dt))[0] + 1

    low = np.append([0], cut)
    high = np.append(cut, len(orig_lc.time))

    smooth_signal = np.zeros(orig_lc.time.size)

    # Smooth each segment separately to avoid artifacts
    for l, h in zip(low, high):
        segment_length = h - l
        # If the segment is too short, just take the median
        if np.any([default_gauss_x > segment_length, segment_length < break_tolerance]):
            smooth_signal[l:h] = np.nanmedian(orig_lc.flux[l:h])
        else:
            # Convolve data
            smooth_signal[l:h] = convolve_teq(orig_lc.flux[l:h], kernel, **kwargs)

    return smooth_signal


def error_estimate(orig_lc: LightCurve, kernel=None, quality=None, boundary=None, **kwargs):
    """
    Estimate the error of a light curve by smoothing the signal, and calculating the
    Mean Absolute Error of the original light curve compared to the smoothed signal.


    :param orig_lc: light curve
    :param kernel: The convolution kernel. Must be 1-Dimensional. Default Gaussian1DKernel
    :param quality: mask array that has usable data for light curve
    :param boundary: flag indicating how to handle boundaries. See convolve or convolve_fft docs.
    :return: return the Median Absolute Error of the data in parts-per million (ppm)
    """

    if quality is None:
        quality = slice(None)

    # Get smooth signal
    smooth_lc = smooth_lightcurve(orig_lc, kernel, boundary=boundary, **kwargs)

    # Calculate residuals, and sigma-clip them
    residuals = orig_lc.flux[quality] - smooth_lc[quality]

    # return the MAE of residuals
    return np.median(np.abs(residuals))


def amplitude_estimate(orig_lc: LightCurve, low=1, high=99, quality=None):
    """
    Estimate the amplitude of a light curve using the low-th, and high-th percentiles of the data.

    :param orig_lc: light curve
    :param high: high percentile value
    :param low: low percentile value
    :return:
    """
    if quality is None:
        quality = slice(None)

    # Amplitude determination; max - min
    flux_min, flux_max = np.percentile(orig_lc.flux[quality], [low, high])
    amplitude = (flux_max - flux_min) / 2
    return amplitude
