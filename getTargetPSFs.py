from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Column
from astropy.stats import sigma_clip
from astroquery.mast import Tesscut
from lightkurve import TessTargetPixelFile, search_targetpixelfile, search_tesscut, TessLightCurve, TessLightCurveFile, MPLSTYLE
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from random import randint, sample
import warnings

percentile=80
radius = 3

def plot_frame(tpf, aperture=None, ax=None, savefn=None, frame=200, show_colorbar=True, **kwargs):
    if not ax:
        ax = plt.subplot(projection=tpf.wcs)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        if aperture=='threshold':
            aperture = tpf.create_threshold_mask()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        tpf.plot(ax=ax, interpolation="nearest", aperture_mask=aperture,cmap="hot", scale='sqrt', frame=frame, show_colorbar=show_colorbar, **kwargs)
    with plt.style.context(MPLSTYLE):
        ax.set_ylabel('Declination')
        ax.set_xlabel('Right Ascension')
    
    # IF want to save
    if savefn:
        plt.gcf().savefig(savefn)

        
if __name__ == '__main__':
    # Import target list
    targets = ascii.read('DataInput/cluster_targets_tic.ecsv')

    sampletpfs = glob('DataOutput/SampleTPFs/*.fits')

    # Make parallel list of target names; convert ID's into integers
    sample_targets = [filename.split('_')[1] for filename in sampletpfs]

    # G Mag Group ranges
    gmin, gmax = targets['G Group'].min(), targets['G Group'].max()

    # Group Apertures
    ap_shape = TessLightCurveFile(sampletpfs[0]).hdu[2].data.shape
    num_groups = gmax-gmin+1
    group_apts = np.zeros((num_groups,) + ap_shape)

    # Find Group Samples, and Calculate Apertures
    for i, group in enumerate(range(gmin, gmax+1)):
        print(f'Collecting Group Samples for GMag={group}')
        grows = targets['G Group'] == group
        group_targets = targets[grows]['TIC ID']

        apertures = []
        for tpf, target in zip(sampletpfs, sample_targets):
            if target in group_targets:
                iapt = TessLightCurveFile(tpf).hdu[2].data
                apertures.append(iapt)


        group_apts[i] = np.sum(apertures, axis=0) > max(6, len(apertures))


    # Loop over targets to plot and save them
    for target in targets:
        ticid = target['TIC ID']
        print(f"Downloading target pixel file for {target['TIC ID']}")
        tpf = search_tesscut(ticid, sector=8).download(cutout_size=8)
        target_apt = group_apts[target['G Group']-gmin]

        print(f"Creating plot for target {target['TIC ID']}")

        # Create percentile aperture
        median_image = np.nanmedian(tpf.flux, axis=0)
        rows, cols = median_image.shape
        x, y = np.ogrid[-rows/2:rows/2, -cols/2:cols/2]
        radius_mask = np.sqrt(x**2 + y**2) < radius
        above_percentile = median_image > np.nanpercentile(median_image, percentile)
        per_mask = above_percentile & radius_mask

        # Create a grid of plots
        nrows, ncols = 2, 2
        dx = 1
        figsize = plt.figaspect(float(dx * nrows) / float(dx * ncols))

        fig, axes = plt.subplots(nrows,ncols, sharex='col', sharey='row', subplot_kw={'projection': tpf.wcs}, 
                                 gridspec_kw={'hspace': 0.2, 'wspace': 0}, figsize=1.5*figsize)
        
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


        # pad = 0.05 # Padding around the edge of the figure
        # xpad, ypad = dx * pad, dx * pad
        # fig.subplots_adjust(left=xpad, right=1-xpad, top=1-ypad, bottom=ypad)

        plt.tight_layout()
        plt.subplots_adjust(top=0.85, left=0.132, bottom=0.1),# wspace=0., hspace=0.2)
        fig.savefig(f"./giantplot/{ticid}_M{target['G Group']}_FFI.jpg")
        plt.close()    

        print(f"Saved plot for target {target['TIC ID']}")

    
#     # Plot and Save FFI w aperture
#     ffiap_savefn = f"./TargetFFI_OldAp/{ticid}_M{target['G Group']}_FFIAP.jpg"
#     plot_frame(tpf, aperture=target_apt, savefn=ffiap_savefn)
    
#     # Plot and Save FFI w/ Threshold Ap    
#     ffiqap_savefn = f"./TargetFFI_ThreshAp/{ticid}_M{target['G Group']}_FFIThAP.jpg"
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore", RuntimeWarning)
#         plot_frame(tpf, aperture='threshold', savefn=ffiqap_savefn)

#     # Plot and Save FFI w/ Percentile Ap
#     ffiperap_savefn = f"./TargetFFI_PercentiledAp/{ticid}_M{target['G Group']}_FFIPerAP.jpg"
#     plot_frame(tpf, aperture=per_mask, savefn=ffiperap_savefn)
    
#     # Plot and Save FFI w/o Aperture
#     ffi_savefn = f"./TargetFFI/{ticid}_M{target['G Group']}_FFI.jpg"
#     plot_frame(tpf, savefn=ffi_savefn)    
    
    
    
