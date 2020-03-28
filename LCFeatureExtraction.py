import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve, convolve_fft
from lightkurve import LightCurve

default_stddev, default_gauss_x = 3, 37


def smooth_lightcurve(orig_lc: LightCurve, kernel=None, break_tolerance=5, convolution='convolve', **kwargs):
    """
    Smooth a light curve with a given kernel through convolution.

    :param orig_lc: light curve to be smoothed.
    :param kernel: The convolution kernel. Must be 1-Dimensional. Default Gaussian1DKernel
    :param break_tolerance: tolerance for which to breake light curve into sub-lightcurves
    :param convolution: convolution function; default astropy.stats..convolution
    :param kwargs: Keyword arguments to pass to the convolution function
    :return:
    """
    if kernel is None:
        # Create kernel
        kernel = Gaussian1DKernel(stddev=default_stddev, x_size=default_gauss_x)
    if convolution == 'convolve':
        convolve_teq = convolve
    elif convolution == 'fft':
        convolve_teq = convolve_fft
    else:
        raise ValueError(f'Wrong convolution technique "{convolution}"')

    # Split the lightcurve into segments by finding large gaps in time
    dt = orig_lc.time[1:] - orig_lc.time[0:-1]
    cut = np.where(dt > break_tolerance * np.nanmedian(dt))[0] + 1

    low = np.append([0], cut)
    high = np.append(cut, len(orig_lc.time))

    smooth_signal = np.zeros(orig_lc.time.size)

    # Smooth each segment separately to avoid artifacts
    for l, h in zip(low, high):
        gap_length = h - l
        # If the segment is too short, just take the median
        if np.any([default_gauss_x > gap_length, gap_length < break_tolerance]):
            smooth_signal[l:h] = np.nanmedian(orig_lc.flux[l:h])
        else:
            # Convolve data
            smooth_signal[l:h] = convolve(orig_lc.flux[l:h], kernel, **kwargs)

    return smooth_signal


def error_estimate(orig_lc: LightCurve, kernel=None, sigma=3,
                   maxiters=3, quality=None, boundary=None, **kwargs):
    """
    Estimate the error of a light curve by smoothing the signal, and calculating the
    Mean Absolute Error of the original light curve compared to the smoothed signal.


    :param orig_lc: light curve
    :param kernel: The convolution kernel. Must be 1-Dimensional. Default Gaussian1DKernel
    :param quality: mask array that has usable data for light curve
    :param boundary: flag indicating how to handle boundaries. See convolve or convolve_fft docs.
    :param sigma: Number of standard deviations to use for both the lower and upper clipping limit.
    :param maxiters: max number of sigma-clipping iterations to perform or None to clip until convergence is achieved
    :return: return the Median Absolute Error of the data in parts-per million (ppm)
    """

    if quality is None:
        quality = slice(None)

    # Get smooth signal
    smooth_lc = smooth_lightcurve(orig_lc, kernel, boundary=boundary, **kwargs)

    # Calculate residuals, and sigma-clip them
    residuals = orig_lc.flux[quality] - smooth_lc[quality]
    # clipped_residuals = sigma_clip(residuals, sigma=sigma, maxiters=maxiters, masked=True)

    # return the MAE of residuals
    return np.median(np.abs(residuals))
    # return np.median(np.abs(residuals)) * 1e6
    # return np.mean(np.abs(residuals)) * 1e6
    # return np.mean(np.abs(clipped_residuals)) * 1e6


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
