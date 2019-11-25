import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits, ascii
from astropy.stats import sigma_clipped_stats
from lightkurve import MPLSTYLE
from lightkurve.utils import plot_image
from matplotlib import patches

from usefulFuncs import ifind_tpfs, getPercentileAp


def find_mag(targets, ticid):
    """
    Simple function to return the magnitude group of a target
    :param targets: master table with targets
    :param ticid: TIC ID of target
    :return: mag group
    """
    return targets['G Group'][targets['TIC ID'] == ticid][0]


percentile = 80
radius = 3
sector = 8
cutsize = 8

titles = 'Pipeline Ap', 'Percentile Ap', 'Threshold Ap'
mask_color = 'pink'

# Get all target pixel files
all_tpfs = ifind_tpfs()

# Import target list
targets = ascii.read('DataInput/cluster_targets_tic.ecsv')

# Get quality mask
with fits.open('DataInput/ClusterQuality/Sector8_Sample.fits.gz') as quality_sample:
    # For Chelsea's quality flags:
    #   0 means good, 1 means bad
    quality_flags = quality_sample[1].data['quality']
    bool_flags = ~quality_flags.astype(bool)

    colors = ['k' if flag else 'red' for flag in bool_flags]

for tpf in all_tpfs:
    # Get percentile aperture before disable default quality flags
    print(f'Calculating apertures for {tpf.targetid}')
    gmag = find_mag(targets, tpf.targetid)

    percAp, above_percentile = getPercentileAp(tpf, both=True, iterative=True)
    threshAp = tpf.create_threshold_mask()

    origPix = tpf.pipeline_mask.sum()
    percPix = percAp.sum()
    threshPix = threshAp.sum()

    # In order to use my own quality flags I must disable current flags
    tpf.quality_mask = np.ones(tpf.quality_mask.shape, dtype=bool)
    tpf.hdu[1].data['QUALITY'] = np.zeros(tpf.hdu[1].data['QUALITY'].shape, dtype=np.int32)

    # Get full lightcurves of each
    lc_orig = tpf.extract_aperture_photometry(aperture_mask='pipeline')
    lc_thresh = tpf.extract_aperture_photometry(aperture_mask=threshAp)
    lc_perc = tpf.extract_aperture_photometry(aperture_mask=percAp)

    # Create copy of light curve to save background-subtracted ones
    lc_origbs = lc_orig.copy()
    lc_threshbs = lc_thresh.copy()
    lc_percbs = lc_perc.copy()

    # There aren't a lot of ways to get a decent background
    # I think a good starting point is percentile aperture
    mean_orig, median_orig, std_orig = sigma_clipped_stats(tpf.flux, sigma=3.0,
                                                           mask=np.broadcast_to(tpf.pipeline_mask, tpf.shape),
                                                           axis=(1, 2))

    mean_thresh, median_thresh, std_thresh = sigma_clipped_stats(tpf.flux, sigma=3.0,
                                                                 mask=np.broadcast_to(threshAp, tpf.shape),
                                                                 axis=(1, 2))

    mean_perc, median_perc, std_perc = sigma_clipped_stats(tpf.flux, sigma=3.0,
                                                           mask=np.broadcast_to(percAp, tpf.shape),
                                                           axis=(1, 2))

    # Background subtract and propagate error (ONLY USING MEAN)
    lc_origbs.flux -= origPix * mean_orig
    lc_origbs.flux_err = np.sqrt(lc_origbs.flux_err ** 2 + (std_orig * origPix) ** 2)

    lc_threshbs.flux -= threshPix * mean_thresh
    lc_threshbs.flux_err = np.sqrt(lc_threshbs.flux_err ** 2 + (std_thresh * threshPix) ** 2)

    lc_percbs.flux -= percPix * mean_perc
    lc_percbs.flux_err = np.sqrt(lc_percbs.flux_err ** 2 + (std_perc * percPix) ** 2)

    aps = [tpf.pipeline_mask, percAp, threshAp]

    frame = 300
    img_extent = (tpf.column, tpf.column + tpf.shape[2],
                  tpf.row, tpf.row + tpf.shape[1])

    with plt.style.context(MPLSTYLE):
        print(f'Creating plot for target {tpf.targetid}')
        fig = plt.figure(figsize=(12, 12))
        grid = plt.GridSpec(4, 3, hspace=0.6, wspace=0.1, bottom=0.05, left=0.08, right=1.,
                            width_ratios=[1, 1, 1], height_ratios=[1, 1, 1, 1])
        top_shared = dict(aspect='equal', projection=tpf.wcs)
        ax1 = fig.add_subplot(grid[0, 0], **top_shared)
        ax2 = fig.add_subplot(grid[0, 1], sharey=ax1, sharex=ax1, **top_shared)
        ax3 = fig.add_subplot(grid[0, 2], sharey=ax1, sharex=ax1, **top_shared)
        ax4 = fig.add_subplot(grid[1, :])
        ax5 = fig.add_subplot(grid[2, :], sharey=ax4, sharex=ax4)
        ax6 = fig.add_subplot(grid[3, :], sharey=ax4, sharex=ax4)
        axes = ax1, ax2, ax3

        ax1.coords[1].set_axislabel('')
        for ax, title, ap in zip(axes, titles, aps):
            plot_image(tpf.flux[frame], ax=ax, title=title, show_colorbar=False, extent=img_extent)
            ax.coords.grid(True, color='white', ls='solid')
            ax.coords[0].set_axislabel('Right Ascension')
            ax.set_ylabel('')
            aperture_mask = tpf._parse_aperture_mask(ap)
            ax.tick_params(axis='x', which='major', labelsize=10)
            for i in range(tpf.shape[1]):
                for j in range(tpf.shape[2]):
                    if aperture_mask[i, j]:
                        ax.add_patch(patches.Rectangle((j + tpf.column, i + tpf.row),
                                                       1, 1, color=mask_color, fill=True,
                                                       alpha=.6))
        this, mask = lc_origbs.remove_outliers(sigma=2.5, return_mask=True)
        this.scatter(ax=ax4, normalize=True, c=[c for c, b in [*zip(colors, ~mask)] if b],
                     show_colorbar=False)
        ax4.set_title('Pipeline Aperture LC')
        this2, mask = lc_percbs.remove_outliers(sigma=2.5, return_mask=True)
        this2.scatter(ax=ax5, normalize=True, c=[c for c, b in [*zip(colors, ~mask)] if b],
                      show_colorbar=False)
        ax5.set_title('Percentile Aperture LC')
        this3, mask = lc_threshbs.remove_outliers(sigma=2.5, return_mask=True)
        this3.scatter(ax=ax6, normalize=True, c=[c for c, b in [*zip(colors, ~mask)] if b],
                      show_colorbar=False)
        ax6.set_title('Threshold Aperture LC')

        ax1.coords[1].set_axislabel('Declination')
        ax4.set_xlabel('')
        ax5.set_xlabel('')

    print(f'Saving plot to LightCurvesPlots/{tpf.targetid}_M{gmag}_LCs.pdf')
    fig.savefig(f'LightCurvesPlots/{tpf.targetid}_M{gmag}_LCs.pdf', dpi=150)
    plt.close()

    # Save LC FITs
    lc_origbs.to_fits(f'./LightCurvesFITS/{tpf.targetid}OR-targ.fits', overwrite=True)
    lc_threshbs.to_fits(f'./LightCurvesFITS/{tpf.targetid}TH-targ.fits', overwrite=True)
    lc_percbs.to_fits(f'./LightCurvesFITS/{tpf.targetid}PC-targ.fits', overwrite=True)
