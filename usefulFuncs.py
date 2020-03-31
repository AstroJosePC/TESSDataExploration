import re
import warnings
from glob import glob
from os import listdir
from os.path import isfile, join

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from lightkurve import MPLSTYLE, TessLightCurveFile, TessTargetPixelFile
from lightkurve import open as open_tpf
from scipy.ndimage import label

reg = re.compile(r'(\d{8,9})([OPT][RCH])')


def extract_lightcurve(tpf: TessTargetPixelFile, aperture_mask, background_mask=None, sigma=3.0,
                       remove_outlier=True, background_sub=False, return_mask=False):
    """
    Extract Background Subtracted Light Curve
    Remove outliers (Optional)

    :param tpf: TargetPixelFile containing all the target data
    :param aperture_mask: A boolean array describing the aperture; `True` pixels are those inside aperture
    :param background_mask: A boolean array describing the bkg aperture: `True` pixels are background pixels
    :param sigma: The number of standard deviations to use for clipping limits
    :param remove_outlier: Whether to remove outliers in returned Light Curve
    :param background_sub: Whether to background subtract the light curve
    :param return_mask: Whether to return a mask indicating which data points were removed
    :return: TessLightCurve object with the processed data
    """

    lc_final = tpf.extract_aperture_photometry(aperture_mask=aperture_mask)

    if background_sub:
        if background_mask is None:
            # If no background mask was given, at least mask the given source
            background_mask = aperture_mask
        else:
            background_mask = ~background_mask

        pix_num = aperture_mask.sum()
        # Get sigma clipped statistics per timestamp for each aperture type; masks source target
        mean_orig, median_orig, std_orig = sigma_clipped_stats(tpf.flux, sigma=sigma,
                                                               mask=np.broadcast_to(background_mask, tpf.shape),
                                                               axis=(1, 2), maxiters=3)

        # Background subtract and propagate error (ONLY USING MEAN)
        lc_final.flux -= pix_num * mean_orig
        lc_final.flux_err = np.sqrt(lc_final.flux_err ** 2 + (std_orig * pix_num) ** 2)

    if remove_outlier:
        lc_final, mask = lc_final.remove_outliers(sigma=sigma, return_mask=True)
        if return_mask:
            return lc_final, mask
    return lc_final


def decode_filename(filepath):
    """
    Function to decode LightCurveFITS filenames:
    {TICID}{Type}.fits

    Types: Threshold (TH), Pipeline (OR), Percentile (PER)

    :param filepath: filepath to FITS file
    :return: ticid, aperture type
    """
    ticid, code = reg.search(filepath).groups()
    if code.upper() == 'TH':
        ap_type = 'Threshold'
    elif code.upper() == 'OR':
        ap_type = 'Pipeline'
    else:
        ap_type = 'Percentile'
    return ticid, ap_type


def getOrigAps(target_table='DataInput/cluster_targets_tic.ecsv',
               orig_samples='DataOutput/SampleTPFs/*.fits'):
    """
    Get original apertures that were manually selected during summer exploration 2019
    :param target_table:
    :param orig_samples:
    :return: group apertures in uint8 array
    """
    # Import target list
    if isinstance(target_table, str):
        targets: Table = ascii.read(target_table)
    else:
        # Assume it is a table
        targets: Table = target_table

    sampletpfs = glob(orig_samples)

    # Make parallel list of target names; convert ID's into integers
    sample_targets = [filename.split('_')[1] for filename in sampletpfs]

    # G Mag Group ranges
    gmin, gmax = targets['G Group'].min(), targets['G Group'].max()

    # Group Apertures
    ap_shape = TessLightCurveFile(sampletpfs[0]).hdu[2].data.shape
    num_groups = gmax - gmin + 1
    group_apts = np.zeros((num_groups,) + ap_shape, dtype=np.uint8)

    # Find Group Samples, and Calculate Apertures
    for i, group in enumerate(range(gmin, gmax + 1)):
        print(f'Collecting Group Samples for GMag={group}')
        grows = targets['G Group'] == group
        group_targets = targets[grows]['TIC ID']

        apertures = []
        for tpf, target in zip(sampletpfs, sample_targets):
            if target in group_targets:
                iapt = TessLightCurveFile(tpf).hdu[2].data
                apertures.append(iapt)

        group_apts[i] = np.sum(apertures, axis=0) > max(6, len(apertures))

    return group_apts


def getPercentileAp(tpf, radius=3, percentile=80, both=False, iterative=False):
    """
    Calculate percentile aperture of source in Target Pixel File
    Very persistent on getting pixels at radius
    :param tpf: target pixel file, preferably small cutout
    :param radius: constraint source to be this distance from the center
    :param percentile: pixel values above this percentile count as source, decreases if iterative
    :param both: Whether to return the final aperture and the only-above-percentile aperture
    :param iterative: whether to perform iterative process if mask has no pixels
    :return: aperture mask
    """
    # Create percentile aperture
    median_image = np.nanmedian(tpf.flux, axis=0)
    rows, cols = median_image.shape
    x, y = np.ogrid[-rows / 2:rows / 2, -cols / 2:cols / 2]
    radius_mask = np.sqrt(x ** 2 + y ** 2) < radius
    per_mask = None

    origPer = percentile
    atLastOnce = False
    minPerc = 5 * origPer / 8
    anyPerc = False

    above_percentile = (median_image > np.nanpercentile(median_image, percentile)) & radius_mask

    if not np.any(above_percentile):
        above_percentile = median_image > np.nanpercentile(median_image, percentile)

    # while (percentile > minPerc) and (iterative or not atLastOnce) and not anyPerc:
    #     above_percentile = median_image > np.nanpercentile(median_image, percentile)
    #     per_mask = above_percentile & radius_mask
    #     anyPerc = np.any(per_mask)
    #     percentile -= 5
    #     atLastOnce = True

    # Return only the contiguous region closest to `region`.
    reference_pixel = (tpf.shape[2] / 2, tpf.shape[1] / 2)
    # First, label all the regions:
    labels = label(above_percentile)[0]
    # For all pixels above threshold, compute distance to reference pixel:
    label_args = np.argwhere(labels > 0)
    distances = [np.hypot(crd[0], crd[1])
                 for crd in label_args - np.array([reference_pixel[1], reference_pixel[0]])]
    # Which label corresponds to the closest pixel?
    closest_arg = label_args[np.argmin(distances)]
    closest_label = labels[closest_arg[0], closest_arg[1]]
    percentile_ap = labels == closest_label

    if both:
        return percentile_ap, above_percentile
    else:
        return percentile_ap


def plot_frame(tpf: TessTargetPixelFile, aperture=None, ax=None, savefn=None, frame=200, show_colorbar=True, **kwargs):
    """
    Function to plot TPF frame
    (optional) If aperture is 'threshold' string, a mask will be generated from
    the TPF native create_threshold_mask() function.

    :param tpf: target pixel file
    :param aperture: array-like, 'pipeline', 'all', 'threshold', or None
        A boolean array describing the aperture such that `True` means
        that the pixel will be used.
        If None or 'all' are passed, all pixels will be used.
        If 'pipeline' is passed, the mask suggested by the official pipeline
        will be returned.
        If 'threshold' is passed, all pixels brighter than 3-sigma above
        the median flux will be used.
    :param ax: if given, the frame will plot on this axis
    :param savefn: safe filename if want to save file
    :param frame: frame # of TPF which to plot
    :param show_colorbar: whether to show color bar
    :param kwargs: keyword arguments to pass to tpf.plot function
    """
    if not ax:
        ax = plt.subplot(projection=tpf.wcs)
    # Set default plotting args
    kwargs['interpolation'] = 'nearest'
    kwargs['cmap'] = 'hot'
    kwargs['scale'] = 'sqrt'

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        tpf.plot(ax=ax, frame=frame, show_colorbar=show_colorbar, aperture_mask=aperture, **kwargs)
    with plt.style.context(MPLSTYLE):
        ax.coords[0].set_axislabel('Right Ascension')
        ax.coords[1].set_axislabel('Declination')
    # IF want to save
    if savefn:
        plt.gcf().savefig(savefn)
    return ax


def find_tpf(ticid, sector=8, cutsize=8):
    fname = f'./TESSCuts/TPF{ticid}_S{sector}C{cutsize}.fits'
    if isfile(fname):
        return open_tpf(fname)
    else:
        print(f'Target Pixel File w/ TICID={ticid} not found')


def find_tpfs(ticids=None, sectors=None, cutsizes=None):
    """
    Retrieve Target Pixel Files from saved cutouts
    :param ticids: list of TIC IDs of targets, string or list
    :param sectors: sectors for which to retrieve, int or list
    :param cutsizes: cut out sizes for which to retrieve, int or list
    :return: list of target pixel file objects
    """
    return [*ifind_tpfs(ticids=ticids, sectors=sectors, cutsizes=cutsizes)]


def ifind_tpfs(ticids=None, sectors=None, cutsizes=None):
    """
    Retrieve Target Pixel Files from saved cutouts; as a generator
    :param ticids: list of TIC IDs of targets, string or list
    :param sectors: sectors for which to retrieve, int or list
    :param cutsizes: cut out sizes for which to retrieve, int or list
    :return: generator of target pixel file objects
    """
    template_fn = r'^TPF{ticid}_S{sec}C{cut}.fits$'
    cuts = secs = r'\d{1,2}'

    if sectors:
        sectorsarr = np.atleast_1d(sectors)
        secs = fr'({"|".join(sectorsarr)})'
    if cutsizes:
        cutsarr = np.atleast_1d(cutsizes)
        cuts = fr'({"|".join(cutsarr)})'

    all_tpfs = listdir('TESSCuts')

    if not ticids:
        ticidsre = r'\d{8,9}'
    else:
        ticidsarr = np.atleast_1d(ticids)
        ticidsre = fr'({"|".join(ticidsarr)})'

    regexp = re.compile(template_fn.format(ticid=ticidsre, sec=secs, cut=cuts))
    for tpf in all_tpfs:
        if regexp.match(tpf):
            tpf_fn = join('.', 'TESSCuts', tpf)
            yield open_tpf(tpf_fn)
