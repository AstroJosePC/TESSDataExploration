"""
This script will take ALL the TPFs and apply the Handpicked, Threshold, and Percentile
apertures to produce light curves, and save them in both FITS and PDF (image).
Producing the light curves and saving them to FITS files is easy, but most of the
code here is creating the perfect plot for showcasing the effets of 
each aperture method on its light curve.

The LightCurve FITS files will saved to a `LightCurvesFITS` folder, and the plots
will be saved to a `LightCurvesPlots` folder.

:param targets_file: Filepath to target list, assumes astropy compatible table format
:param ticid_col: Column in input table that has the targets TIC ID
:param mag_col: Column string that contains the grouping magnitude to assign apertures
:param quality_path: Path to Quality flags FITS file, in my case it followed Chelsea's quality flags.

"""
from os import mkdir
from os.path import isdir

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits, ascii
from lightkurve import MPLSTYLE
from lightkurve.utils import plot_image
from matplotlib import patches
from matplotlib.lines import Line2D

from usefulFuncs import ifind_tpfs, getPercentileAp, extract_lightcurve

targets_file = 'DataInput/cluster_targets_tic.ecsv'
tic_col = 'TIC ID'
mag_col = 'G Group'
quality_path = 'DataInput/ClusterQuality/Sector8_Sample.fits.gz'
percentile = 80
radius = 3
sector = 8
cutsize = 8
sigma = 3.0
mask_color = 'pink'


def find_mag(targets, ticid):
    """
    Simple function to return the magnitude group of a target
    :param targets: master table with targets
    :param ticid: TIC ID of target
    :return: mag group
    """
    return targets[mag_col][targets[tic_col] == ticid][0]


# Get all target pixel files
all_tpfs = ifind_tpfs()

# Import target list
targets = ascii.read(targets_file)

# Create folder to save LightCurve FITS
if not isdir('LightCurvesFITS'):
    mkdir('LightCurvesFITS')

# Create folder to save LightCurve plots
if not isdir('LightCurvesPlots'):
    mkdir('LightCurvesPlots')

# Get quality mask
with fits.open(quality_path) as quality_sample:
    # For Chelsea's quality flags:
    #   0 means good, 1 means bad
    quality_flags = quality_sample[1].data['quality']
    bool_flags = ~quality_flags.astype(bool)

    colors = ['k' if flag else 'red' for flag in bool_flags]

lc_legend_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor='k', label='Good Quality', markersize=7),
                     Line2D([0], [0], marker='o', color='w', markerfacecolor='red', label='Bad Quality', markersize=7)]

for tpf in all_tpfs:
    # Get percentile aperture before disable default quality flags
    print(f'Calculating apertures for {tpf.targetid}')

    # Get magnitude of star (G-band)
    gmag = find_mag(targets, tpf.targetid)

    # get the percentile, and threshold apertures
    percAp, above_percentile = getPercentileAp(tpf, both=True, iterative=True, percentile=percentile)
    threshAp = tpf.create_threshold_mask()

    origPix = tpf.pipeline_mask.sum()
    percPix = percAp.sum()
    threshPix = threshAp.sum()

    # In order to use my own quality flags I must disable current flags
    # Re-define the quality mask, and rewrite quality array in data fits
    tpf.quality_mask = np.ones(tpf.quality_mask.shape, dtype=bool)
    tpf.hdu[1].data['QUALITY'] = np.zeros(tpf.hdu[1].data['QUALITY'].shape, dtype=np.int32)

    # Extract sigma-clipped background subtracted light curves
    lc_origbs = extract_lightcurve(tpf, tpf.pipeline_mask, background_sub=True, remove_outlier=False)
    lc_threshbs = extract_lightcurve(tpf, tpf.threshAp, background_sub=True, remove_outlier=False)
    lc_percbs = extract_lightcurve(tpf, tpf.percAp, background_sub=True, remove_outlier=False)

    # few parameters for plot_image function; for plotting the frame cutout
    frame = 300
    img_extent = (tpf.column, tpf.column + tpf.shape[2],
                  tpf.row, tpf.row + tpf.shape[1])

    with plt.style.context(MPLSTYLE):
        # For the plot creation we use a GridSpec in order to 'facilitate' the organizing
        # The top three axes (ax1, ax2, ax3) are the target frame cutout's aperture showcase
        # The three long plots are the respective light curves in the same order.
        print(f'Creating plot for target {tpf.targetid}')
        fig = plt.figure(figsize=(12, 12))
        grid = plt.GridSpec(4, 3, hspace=0.58, wspace=0.1, bottom=0.05, left=0.065, right=1.,
                            top=0.99, width_ratios=[1, 1, 1], height_ratios=[1, 1, 1, 1])
        top_shared = dict(aspect='equal', projection=tpf.wcs)
        ax1 = fig.add_subplot(grid[0, 0], **top_shared)
        ax2 = fig.add_subplot(grid[0, 1], sharey=ax1, sharex=ax1, **top_shared)
        ax3 = fig.add_subplot(grid[0, 2], sharey=ax1, sharex=ax1, **top_shared)

        ax5 = fig.add_subplot(grid[2, :])  # sharey=ax4, sharex=ax4)
        ax4 = fig.add_subplot(grid[1, :], sharey=ax5)  # sharex=ax5)
        ax6 = fig.add_subplot(grid[3, :], sharey=ax5)  # sharex=ax5)
        axes = ax1, ax2, ax3

        ax1.coords[1].set_axislabel('')

        # Iterate over aperture photometry methods, and plot frames
        titles = 'Handpicked Ap', 'Percentile Ap', 'Threshold Ap'
        aps = [tpf.pipeline_mask, percAp, threshAp]
        for ax, title, ap in zip(axes, titles, aps):
            plot_image(tpf.flux[frame], ax=ax, title=title, show_colorbar=False, extent=img_extent)
            ax.coords.grid(True, color='white', ls='solid')
            ax.coords[0].set_axislabel('Right Ascension')
            ax.set_ylabel('')
            aperture_mask = tpf._parse_aperture_mask(ap)
            ax.tick_params(axis='both', which='major', labelsize=10, color='w')
            for i in range(tpf.shape[1]):
                for j in range(tpf.shape[2]):
                    if aperture_mask[i, j]:
                        ax.add_patch(patches.Rectangle((j + tpf.column, i + tpf.row),
                                                       1, 1, color=mask_color, fill=True,
                                                       alpha=.6))

        # Plot each light curve; remove outliers, create good/bad quality-color mask, 
        # and finally make scatter plots for each axes
        lc_orig_clean, mask_orig = lc_origbs.remove_outliers(sigma=sigma, return_mask=True)
        color_mask = [c for c, b in [*zip(colors, ~mask_orig)] if b]
        lc_orig_clean.scatter(ax=ax4, normalize=True, c=color_mask, show_colorbar=False)
        ax4.set_title('Handpicked Aperture LC')

        lc_perc_clean, mask_perc = lc_percbs.remove_outliers(sigma=sigma, return_mask=True)
        color_mask = [c for c, b in [*zip(colors, ~mask_perc)] if b]
        lc_perc_clean.scatter(ax=ax5, normalize=True, c=color_mask, show_colorbar=False)
        ax5.set_title('Percentile Aperture LC')

        lc_thresh_clean, mask_thresh = lc_threshbs.remove_outliers(sigma=sigma, return_mask=True)
        color_mask = [c for c, b in [*zip(colors, ~mask_thresh)] if b]
        lc_thresh_clean.scatter(ax=ax6, normalize=True, c=color_mask, show_colorbar=False)
        ax6.set_title('Threshold Aperture LC')

        # Setting Y Limits; further clip data points to reduce Y-axis window
        # This method is very crude and not very logical, I only employed it cause it worked.
        zthis = lc_orig_clean.remove_outliers(sigma=2.5).flux
        zthis2 = lc_perc_clean.remove_outliers(sigma=2.5).flux
        zthis3 = lc_thresh_clean.remove_outliers(sigma=2.5).flux

        zthis /= zthis.mean()
        zthis2 /= zthis2.mean()
        zthis3 /= zthis3.mean()

        ylim_min = np.mean([zthis.min(), zthis2.min(), zthis3.min()])
        ylim_max = np.mean([zthis.max(), zthis2.max(), zthis3.max()])

        ax5.set_ylim(ylim_min, ylim_max)

        ax1.coords[1].set_axislabel('Declination')
        ax4.set_xlabel('')
        ax5.set_xlabel('')

    print(f'Saving plot to LightCurvesPlots/{tpf.targetid}_M{gmag}_LCs.pdf')

    ax4.legend(handles=lc_legend_handles, loc='upper right',
               bbox_to_anchor=(1.0, 1.23), ncol=2, fancybox=True, shadow=False)
    # plt.show()
    fig.savefig(f'LightCurvesPlots/{tpf.targetid}_M{gmag}_LCs.pdf', dpi=300)
    plt.close()

    # Save LC FITs
    lc_origbs.to_fits(f'./LightCurvesFITS/{tpf.targetid}OR-targ.fits', overwrite=True)
    lc_threshbs.to_fits(f'./LightCurvesFITS/{tpf.targetid}TH-targ.fits', overwrite=True)
    lc_percbs.to_fits(f'./LightCurvesFITS/{tpf.targetid}PC-targ.fits', overwrite=True)
