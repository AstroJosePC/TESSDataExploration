from glob import iglob, glob
from os.path import join, isdir, isfile
from os import mkdir
from tqdm import tqdm
from warnings import catch_warnings, simplefilter

import astropy.units as u
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from astropy.table import Table, QTable
from astropy.io.fits.verify import VerifyWarning
from lightkurve import MPLSTYLE
from lightkurve import open as open_lc
from scipy.signal import argrelmax
import numpy as np

from k2spin import prot
from usefulFuncs import decode_filename
from LCFeatureExtraction import error_estimate


def find_mag(targets, ticid):
    """
    Simple function to return the magnitude group of a target
    :param targets: master table with targets
    :param ticid: TIC ID of target
    :return: mag group
    """
    return targets['G Group'][targets['TIC ID'] == ticid][0]


src_lcfs = 'LightCurvesFITS/*.fits'
# fits_paths = iglob(src_lcfs)
fits_paths = glob(src_lcfs)

# Import target list
targets = ascii.read('DataInput/cluster_targets_tic.ecsv')
ticids = targets['TIC ID']

# Create directories for outputs
pgs_savepath = 'DraftPeriodograms'
pgs_fits_savepath = 'PGs_FITS'

if not isdir(pgs_savepath):
    mkdir(pgs_savepath)
if not isdir(pgs_fits_savepath):
    mkdir(pgs_fits_savepath)

# PROGRAM PARAMETERS
outlier_sigma = 2.5

# Parameters for LS & bootstrap algorithm
threshold = 1.
prot_lims = (0.1, 28.)

# Structure: fund_periods[`TICID`][`AP_TYPE`]; Initialize
# Types: Threshold (TH), Pipeline (OR), Percentile (PER)
aperture_types = 'Threshold', 'Pipeline', 'Percentile'
lc_paths = []

# force_targets = ['93269120', '93014257', '92583560', '93549309']
force_targets = []

# Get quality mask
with fits.open('DataInput/ClusterQuality/Sector8_Sample.fits.gz') as quality_sample:
    # For Chelsea's quality flags:
    #   0 means good, 1 means bad
    quality_flags = quality_sample[1].data['quality']
    good_mask = ~quality_flags.astype(bool)

# Iterave over all targets, and find its light curve with lowest noise;
# Also, filter out those we do not need to calculate
for ticid in ticids:
    # By sorting them, I get: Handpicked Ap, Percentile Ap, Threshold Ap
    apt_paths = sorted([fits_path for fits_path in fits_paths if ticid in fits_path])
    sigmas = np.zeros(len(apt_paths))

    for i, apt_path in enumerate(apt_paths):
        clipped_lc, clipped_mask = open_lc(apt_path).get_lightcurve('FLUX').remove_outliers(sigma=outlier_sigma,
                                                                                            return_mask=True)
        sigmas[i] = error_estimate(clipped_lc, quality=good_mask[~clipped_mask])
    
    # Get the index of the path to light curve with lowest noise
    idx = sigmas.argmin()
    
    # Decode filename
    ap_type = decode_filename(apt_paths[idx])[1]
    gmag = find_mag(targets, ticid)
    
    # Generate save paths
    pg_savepath = join(pgs_savepath, f'{ticid}PG-{ap_type}-M{gmag}.pdf')
    pg_fits_savepath = f'{pgs_fits_savepath}/{ticid}PG-{ap_type}-M{gmag}.fits'
    
    # If periodogram info exists, and we aren't forcing the target, then skip file
    if isfile(pg_savepath) and isfile(pg_fits_savepath) and (ticid not in force_targets):
        continue
    else:
        lc_paths.append(apt_paths[idx])


# Setup progress bar as context manager
with tqdm(total=len(lc_paths), unit='pg') as pbar:
    # Iterate over light curve paths, calculate their periodograms, and save data
    for i, lc_path in enumerate(lc_paths):
        ticid, ap_type = decode_filename(lc_path)
        gmag = find_mag(targets, ticid)
        
        # Pre-generate save paths
        pg_savepath = join(pgs_savepath, f'{ticid}PG-{ap_type}-M{gmag}.pdf')
        pg_fits_savepath = f'{pgs_fits_savepath}/{ticid}PG-{ap_type}-M{gmag}.fits'
        
        # Get J-K temp proxy
        j_mag = targets['J'][targets['TIC ID'] == ticid][0]
        k_mag = targets['K'][targets['TIC ID'] == ticid][0]
        j_k = j_mag - k_mag

        # Update progress bar, import Light Curve, and calculate periodogram
        pbar.set_description(f'{ticid}-{ap_type}')
        lc = open_lc(lc_path).get_lightcurve('FLUX').remove_outliers(sigma=outlier_sigma)
        
        # Run LombScarlge Periodogram
        pbar.set_postfix(Status='Running LS')
        fund_period, fund_power, periods_to_test, periodogram, aliases, sigmas = prot.run_ls(lc.time, lc.flux, lc.flux_err,
                                                                                             threshold, prot_lims=prot_lims,
                                                                                             run_bootstrap=True)
        
        pbar.set_postfix(Status='Saving PG FITS')
        meta_info = {'TICID':ticid, 'aperture':ap_type, 'gmag':gmag, 'jkmag':j_k, 'sig99p':sigmas[0]}
        with catch_warnings():
            simplefilter('ignore', VerifyWarning)
            pg_table = QTable(data=(periods_to_test * u.day, periodogram), 
                              names=('period', 'power'), meta=meta_info)
            pg_table.write(pg_fits_savepath, overwrite=True)
        pbar.set_postfix(Status=f'Found period={fund_period:.2f} days')

        # Find local maxima in data; most credible periods
        peaks = argrelmax(periodogram)[0]

        # Remove any peak below threshold
        acceptable_peaks = periodogram[peaks] > sigmas[0]
        new_peaks = peaks[acceptable_peaks]
        period_peaks = periods_to_test[new_peaks]
        
        pbar.set_postfix(Status='Making plot')
        with plt.style.context(MPLSTYLE):
            fig, ax = plt.subplots(figsize=(12, 8))

            pg_curve = ax.plot(periods_to_test, periodogram, 'b-', lw=1.5)
            ax.set_xlabel('Periods [days]', fontsize=18)
            ax.set_ylabel('Normalized Power', fontsize=18)
            ax.set_title(f'Periodogram for TICID {ticid}', fontsize=20)
            threshold_line = ax.axhline(sigmas[0], color='green', linestyle='-.', lw=1.5)

            for period_peak in period_peaks:
                red_dash = ax.axvline(period_peak, ls='--', c='#eb6060', lw=1.5)  # light red shade

            main_red = ax.axvline(fund_period, c='#ad2d2d', lw=1.5)  # dark red shade

            ax.set_xlim(0, 24)
            curves = [pg_curve[0], red_dash, threshold_line, main_red]
            curve_labels = ['Periodogram', 'Peaks', 'Threshold', 'Best Period']

            ax.legend(curves, curve_labels, fontsize=14)
            ax.text(fund_period + 0.5, fund_power, rf'Period={fund_period:0.2f}d', fontsize=14)
            ax.tick_params(axis='both', labelsize=14)

        plt.savefig(pg_savepath, dpi=150)
        plt.close()
        
        pbar.set_postfix(Status=f'Plot saved as {ticid}PG-{ap_type}-M{gmag}.pdf')
        pbar.update()