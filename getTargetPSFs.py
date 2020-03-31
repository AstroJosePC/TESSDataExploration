"""
Create a 2x2 grid of plots showcasing the aperture selection methods for every IC2391 target.

The plots will be saved in TargetFFI folder, with the following name format:  {TIC ID}_M{G-mag}_FFI.jpg
where G-mag is the rounded G-band magnitude of the star.
"""

from os import mkdir
from os.path import isdir

import matplotlib.pyplot as plt
from astropy.io import ascii
from lightkurve import MPLSTYLE

from usefulFuncs import find_tpf, plot_frame, getPercentileAp, getOrigAps

percentile = 80
radius = 3
sector = 8
cutsize = 8

if __name__ == '__main__':
    # Import target list
    targets = ascii.read('DataInput/cluster_targets_tic.ecsv')

    group_apts = getOrigAps(targets).astype(bool)
    gmin = targets['G Group'].min()

    # Create folder to save FFI 
    if not isdir('TargetFFI'):
        mkdir('TargetFFI')

    # Loop over targets to plot and save them
    for ticid, Ggroup in targets[['TIC ID', 'G Group']]:
        print(f"Downloading target pixel file for {ticid}")
        tpf = find_tpf(ticid, sector=sector, cutsize=cutsize)
        if not tpf:
            print('\tSkipping Target Pixel File...')
            continue
        target_apt = group_apts[Ggroup - gmin]

        print(f"Creating plot for target {ticid}")

        # Create percentile aperture
        per_mask = getPercentileAp(tpf, percentile)

        # Create a grid of plots
        nrows, ncols = 2, 2
        dx = 1
        figsize = plt.figaspect(float(dx * nrows) / float(dx * ncols))

        fig, axes = plt.subplots(nrows, ncols, sharex='col', sharey='row', subplot_kw={'projection': tpf.wcs},
                                 gridspec_kw={'hspace': 0.2, 'wspace': 0}, figsize=1.5 * figsize)

        with plt.style.context(MPLSTYLE):
            img_title = 'Full Frame Image\nTarget ID: {}'.format(tpf.targetid, tpf.cadenceno[200])
            plt.suptitle(img_title, fontsize='xx-large', fontweight='bold')

        ((ax1, ax2), (ax3, ax4)) = axes

        plot_frame(tpf, ax=ax1, show_colorbar=False)
        plot_frame(tpf, ax=ax2, aperture=target_apt, show_colorbar=False)
        plot_frame(tpf, ax=ax3, aperture='threshold', show_colorbar=False)
        plot_frame(tpf, ax=ax4, aperture=per_mask, show_colorbar=False)

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

        with plt.style.context(MPLSTYLE):
            ax1.set_title('Frame w/o Ap')
            ax2.set_title('Handpicked Ap')
            ax3.set_title('Threshold Ap')
            ax4.set_title('Percentile Ap')

        plt.tight_layout()
        plt.subplots_adjust(top=0.85, left=0.132, bottom=0.1)
        fig.savefig(f"./TargetFFI/{ticid}_M{Ggroup}_FFI.jpg")
        plt.close()

        print(f"Saved plot for target {ticid}")
