import re
from glob import iglob
from os.path import join, isfile

import matplotlib.pyplot as plt
from astropy.io import ascii
from lightkurve import MPLSTYLE
from lightkurve import open as open_lc
from scipy.signal import argrelmax

from k2spin import prot


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
fits_paths = iglob(src_lcfs)
reg = re.compile(r'(\d{8,9})([OPT][RCH])')

# Import target list
targets = ascii.read('DataInput/cluster_targets_tic.ecsv')

pgs_savepath = 'DraftPeriodograms'

# num_terms is the number of terms for which the fourier expansion is done (how complex are the oscillations?)
num_terms = 1
# Min_period is the minimum detectable period
min_period = 0.1
# max_period is the maximum detectable period
max_period = 28.

ov_sampling = 150
un_sampling = 10

power_cutoff = 75  # percentile

outlier_sigma = 2.5

# Parameters for LS & bootstrap algorithm
threshold = 1.
prot_lims = (0.1, 28.)

fund_periods = dict()
all_ticids = []
periodograms = dict()
ls_results = dict()
J_K = dict()

force_targets = ['93269120', '93014257', '92583560', '93549309']


def decode_filename(filepath):
    ticid, code = reg.search(filepath).groups()
    if code.upper() == 'TH':
        ap_type = 'Threshold'
    elif code.upper() == 'OR':
        ap_type = 'Pipeline'
    else:
        ap_type = 'Percentile'
    return ticid, ap_type


i = 0
for i, fits_path in enumerate(fits_paths):
    # fits_path = fits_paths[25]
    ticid, ap_type = decode_filename(fits_path)
    gmag = find_mag(targets, ticid)

    savepath = join(pgs_savepath, f'{ticid}PG-{ap_type}-M{gmag}.pdf')

    if not force_targets:
        if (gmag > 14) or isfile(savepath):
            continue
    elif ticid not in force_targets:
        continue

    # if (gmag > 14) or isfile(savepath) and (ticid not in force_targets):
    # if gmag > 14:
    #     continue

    j_k = targets['J'][targets['TIC ID'] == ticid][0] - targets['K'][targets['TIC ID'] == ticid][0]

    print(f'Creating periodogram for {ticid}-{ap_type}')
    lcf = open_lc(fits_path)
    lc = lcf.get_lightcurve('FLUX').remove_outliers(sigma=outlier_sigma)

    fund_period, fund_power, periods_to_test, periodogram, aliases, sigmas = prot.run_ls(lc.time, lc.flux, lc.flux_err,
                                                                                         threshold, prot_lims=prot_lims,
                                                                                         run_bootstrap=True)
    print(f'Found best period:\t{fund_period}')
    ikey = f'{ticid}-{ap_type}'

    if ticid not in all_ticids:
        all_ticids.append(ticid)

    J_K[ticid] = j_k
    ls_results[ikey] = fund_period, fund_power, aliases, sigmas
    periodograms[ikey] = periods_to_test, periodogram

    # Find local maxima in data; most credible periods
    peaks = argrelmax(periodogram)[0]

    # Remove any peak below threshold
    acceptable_peaks = periodogram[peaks] > sigmas[0]
    new_peaks = peaks[acceptable_peaks]
    period_peaks = periods_to_test[new_peaks]

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

    print(f'Saving figure to {savepath}')
    plt.savefig(savepath, dpi=150)
    plt.close()
#     break

print('Done')
# if not np.any(np.isclose(fund_period, prot_lims)):
#     if ticid not in ticids:
#         J_K.append(j_k)
#         ticids.append(ticid)
#         fund_periods[ticid] = [fund_period,]
#     else:
#         fund_periods[ticid].append(fund_period)
#
# plt.show()

# plt.scatter(J_K, fund_periods)
# plt.title('found periods in IC2391 <14Mags')
# plt.savefig('periods.jpg', dpi=150)
# plt.close()

# th_periods = []
# pipe_periods = []
# perc_periods = []
# ordered_mags = []
#
# for ticid in all_ticids:
#     th_key = f'{ticid}-Threshold'
#     if th_key in ls_results:
#         th_periods.append(ls_results[th_key][0])
#     else:
#         th_periods.append(np.nan)
#
#     pipe_key = f'{ticid}-Pipeline'
#     if pipe_key in ls_results:
#         pipe_periods.append(ls_results[pipe_key][0])
#     else:
#         pipe_periods.append(np.nan)
#
#     perc_key = f'{ticid}-Percentile'
#     if perc_key in ls_results:
#         perc_periods.append(ls_results[perc_key][0])
#     else:
#         perc_periods.append(np.nan)
#
#     ordered_mags.append(J_K[ticid])


# # Undersample periodogram to get a feeling on main local maxima
# pg = lc.to_periodogram(method='lombscargle', minimum_period=min_period*u.day, maximum_period=max_period*u.day,
#                        nterms=num_terms, oversample_factor=un_sampling)

# # Find local maxima in data; most credible periods
# peaks = argrelmax(pg.power)[0]
#
# top_10_peaks = peaks[pg.power[peaks].argsort()[:-10:-1]]
#
# period_peaks = pg.period[top_10_peaks]
#
# # Plot periodogram in period domain!!
# ax = pg.plot(scale='linear', view='period', title='Periodogram')
#
# for period_peak in period_peaks:
#     plt.axvline(period_peak.value, ls='--', c='red')
#
# max_period_plot = min(period_peaks.max().value + 1., pg.period.max().value)
# max_power_plot = min(pg.power[top_10_peaks].max().value + 20., pg.power.max().value)
# plt.xlim(0, max_period_plot)
# plt.ylim(0, max_power_plot)
# plt.show()
