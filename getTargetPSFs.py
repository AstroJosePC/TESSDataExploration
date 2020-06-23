"""
Create a 2x2 grid of plots showcasing the aperture selection methods for every IC2391 
target; bare, handpicked, threshold, and percentile. This is good to compare the
aperture sizes/shapes with each other.

This script is harder to make portable since it requires a way to extract handpicked 
apertures. Take a look at the notes on how to make handpicked apertures; at least the
way I made them.

The plots will be saved in TargetFFI folder, with the following name format:  {TIC ID}_M{Mag}_FFI.jpg
where Mag is the rounded #-band magnitude of the star used to group the handpicked magnitudes.

There are a few paramters to take care of:

:param targets_file: Filepath to target list, assumes astropy compatible table format
:param percentile: percentile value to use for percentile aperture generation
:param sector: sector number to download data from; default Sector 8.
:param cutsize: frame cutout size to get; default 8x8
:param ticid_col: Column in input table that has the targets TIC ID
:param mag_col: Column string that contains the grouping magnitude to assign apts
"""

from os import mkdir
from os.path import isdir

import matplotlib.pyplot as plt
from astropy.io import ascii
from lightkurve import MPLSTYLE

from usefulFuncs import find_tpf, plot_frame, getPercentileAp, getOrigAps

targets_file = 'DataInput/cluster_targets_tic.ecsv'
percentile = 80
sector = 8
cutsize = 8
ticid_col = 'TIC ID'
mag_col = 'G Group'

if __name__ == '__main__':
    # Import target list
    targets = ascii.read(targets_file)

    # Extract handpicked apertures 
    group_apts = getOrigAps(targets).astype(bool)
    mag_min = targets[mag_col].min()

    # Create folder to save FFI 
    if not isdir('TargetFFI'):
        mkdir('TargetFFI')

    # Loop over targets to plot and save them
    for ticid, mag in targets[[ticid_col, mag_col]]:
        print(f"Downloading target pixel file for {ticid}")
        tpf = find_tpf(ticid, sector=sector, cutsize=cutsize)
        if not tpf:
            print('\tSkipping Target Pixel File...')
            continue

        # Get the target aperture by indexing it with the magnitude orders
        target_apt = group_apts[mag - mag_min]

        print(f"Creating plot for target {ticid}")

        # Create percentile aperture
        per_mask = getPercentileAp(tpf, percentile)

        # Create a grid of plots
        nrows, ncols = 2, 2
        dx = 1
        figsize = plt.figaspect(float(dx * nrows) / float(dx * ncols))
        fig, axes = plt.subplots(nrows, ncols, sharex='col', sharey='row', 
                                 subplot_kw={'projection': tpf.wcs}, 
                                 gridspec_kw={'hspace': 0.2, 'wspace': 0}, 
                                 figsize=1.5 * figsize)

        # Setup a suptitle with lightkurve style
        with plt.style.context(MPLSTYLE):
            img_title = 'Full Frame Image\nTarget ID: {}'.format(tpf.targetid, 
                                                                 tpf.cadenceno[200])
            plt.suptitle(img_title, fontsize='xx-large', fontweight='bold')

        # Separate axes to pass, and plot each cutout frame into them with their apts
        ((ax1, ax2), (ax3, ax4)) = axes
        plot_frame(tpf, ax=ax1, show_colorbar=False)
        plot_frame(tpf, ax=ax2, aperture=target_apt, show_colorbar=False)
        plot_frame(tpf, ax=ax3, aperture='threshold', show_colorbar=False)
        plot_frame(tpf, ax=ax4, aperture=per_mask, show_colorbar=False)

        # Setup RA/Dec Axes for informational purposes
        for ax in axes.flat:
            ra = ax.coords[0]
            dec = ax.coords[1]

            if not ax.is_first_col():
                dec.set_ticks_visible(False)
                dec.set_ticklabel_visible(False)
                dec.set_axislabel('')

            if not ax.is_last_row():
                ra.set_ticks_visible(False)
                ra.set_ticklabel_visible(False)
                ra.set_axislabel('')

        # Setup titles 
        with plt.style.context(MPLSTYLE):
            ax1.set_title('Frame w/o Ap')
            ax2.set_title('Handpicked Ap')
            ax3.set_title('Threshold Ap')
            ax4.set_title('Percentile Ap')

        plt.tight_layout()
        plt.subplots_adjust(top=0.85, left=0.132, bottom=0.1)
        fig.savefig(f"./TargetFFI/{ticid}_M{mag}_FFI.jpg")
        plt.close()

        print(f"Saved plot for target {ticid}")
