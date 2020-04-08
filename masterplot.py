from glob import iglob, glob
from os.path import join, isfile
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u
from astropy.io import ascii, fits
import numpy as np
from lightkurve import MPLSTYLE, LombScarglePeriodogram
from lightkurve import open as open_lc
from scipy.signal import argrelmax
from usefulFuncs import decode_filename
from LCFeatureExtraction import error_estimate

targets = ascii.read('DataInput/cluster_targets_tic.ecsv')
ticids = targets['TIC ID']
solar_jk = (0.34, 0.49)
fast_period = 0.4
med_period = 6
slow_period = 8
mark_size = 150.0

aperture_types = 'Threshold', 'Pipeline', 'Percentile'
src_pgfs = 'PGs_FITS/*.fits'
pgfs_paths = glob(src_pgfs)
num_pfs = len(pgfs_paths)

fund_periods = np.zeros(num_pfs)
fund_powers = np.zeros(num_pfs)
low_sigmas = np.zeros(num_pfs)
all_ticids = [] * num_pfs
all_j_k = np.zeros(num_pfs)
all_aptypes = [] * num_pfs
group_gmags = np.zeros(num_pfs, dtype=int)

for i, pgf_path in enumerate(pgfs_paths):
	with fits.open(pgf_path) as pgf:
		max_pow = pgf[1].data['power'].argmax()
		fund_periods[i] = pgf[1].data['period'][max_pow]
		all_j_k[i] = pgf[1].header['jkmag']
		low_sigmas[i] = pgf[1].header['sig99p']
		group_gmags[i] = pgf[1].header['gmag']
		all_ticids.append(pgf[1].header['TICID'])
		all_aptypes.append(pgf[1].header['aperture'])
    
# This master table contains all the data we need to make master plot
master_table = Table(data=(all_ticids, group_gmags, all_aptypes, all_j_k, fund_periods, fund_powers, low_sigmas),
                     names=('tic id', 'g group' ,'ap type', 'j-k', 'period', 'power', 'sigma'))

master_table.write('master_table.fits', overwrite=True)

indx_points = (0.5 < fund_periods)  & (fund_periods < 28.0)
periods = fund_periods[indx_points]
mags = all_j_k[indx_points]
num_periods = periods.size
max_period = periods.argmax()

# Convert periods to angular velocity in units of solar angular velocity
# omega_sun=2.87 × 10^−6 s^−1. From Gallet&Bouvier 2013 -->  Weber & Davis 1967
days2secs = 86400
omega_sun = 2.87e-6  # 1/s
Period2Omega = 2 * np.pi / (days2secs * omega_sun)
norm_omegas = Period2Omega / periods 

with plt.style.context(MPLSTYLE):
	# Create Log plot
	fig, ax = plt.subplots(figsize=(12, 8))
	ax.scatter(mags, norm_omegas, s=mark_size, marker='+')
	ax.set_xlabel('J - K', fontsize=22)
	ax.set_ylabel(r'$\Omega_* / \Omega_\odot$', fontsize=22)
	ax.set_yscale('log')
	ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x,pos:f'{int(x)}'))
	ax.tick_params(top=True, which='both', labelsize=16)
	ax.tick_params(which='major', length=10, width=2)
	ax.tick_params(which='minor', length=6, width=1)

	# Get autoscaled axes limits, and set them on secondary axis
	norm_omega_lims = np.array(ax.get_ylim())
	period_max, period_min = np.abs(Period2Omega / norm_omega_lims)

	# Get secondary axis
	ax2 = ax.twinx()
	ax2.set_ylabel(r'$\mathrm{P_{rot}}$ (d)', rotation=-90, fontsize=22, labelpad=25)
	ax2.set_yscale('log')
	ax2.set_ylim(period_max, period_min)
	ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x,pos:(f'{x:.1f}') if (x < 1) else (f'{int(x)}')))
	ax2.tick_params(labelsize=16, which='major', length=10, width=2)
	ax2.tick_params(which='minor', length=6, width=1)

	ax.set_title(f'{num_periods} Periods from IC2391', fontsize=26)
	plt.tight_layout()
	plt.savefig('masterplot_log.pdf', dpi=150)
	plt.close()


	# Create regular plot
	fig, ax = plt.subplots(figsize=(12, 8))
	ax.scatter(mags, periods, s=mark_size, marker='+')
	ax.set_xlabel('J - K', fontsize=22)
	ax.set_ylabel(r'$\mathrm{P_{rot}}$ (d)', fontsize=22)
	ax.tick_params(top=True, which='both', labelsize=16)
	ax.tick_params(which='major', length=10, width=2)
	ax.tick_params(which='minor', length=6, width=1)

	ax.set_title(f'{num_periods} Periods from IC2391', fontsize=26)
	plt.tight_layout()
	plt.savefig('masterplot.pdf', dpi=150)
	plt.close()

solar_like = (0.24 < mags) & (mags < 0.65)
solar_periods = periods[solar_like]
solar_mags = mags[solar_like]
norm_omegas_solar = Period2Omega / solar_periods
num_solar = solar_periods.size

with plt.style.context(MPLSTYLE):
	plot_title = f'{num_solar} Solar Periods from IC2391 ($0.75< M_*/M_\odot <1.25$)'
	plot_title += f'\nFastest Rotator: {solar_periods.min():.2f} d OR {norm_omegas_solar.max():.1f} $\Omega_\odot$'

	# Create Log plot
	fig, ax = plt.subplots(figsize=(12, 8))
	ax.scatter(solar_mags, norm_omegas_solar, s=mark_size, marker='+')
	ax.axhline(Period2Omega/fast_period, lw=1.0, c='blue', label='ZAMS Fast Rot')
	# ax.axhline(Period2Omega/med_period, lw=1.0, c='green', label='ZAMS Med Rot')
	ax.axhline(Period2Omega/slow_period, lw=1.0, c='red', label='ZAMS Slow Rot')
	ax.axvline(solar_jk[0], ls='--', lw=1.0, c='#a15f43')
	ax.axvline(solar_jk[1], ls='--', lw=1.0, c='#a15f43', label=r'0.9-1.1 $M_\odot$')

	ax.set_xlabel('J - K', fontsize=22)
	ax.set_ylabel(r'$\Omega_* / \Omega_\odot$', fontsize=22)
	ax.set_yscale('log')
	ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x,pos:f'{int(x)}'))
	ax.tick_params(top=True, which='both', labelsize=16)
	ax.tick_params(which='major', length=10, width=2)
	ax.tick_params(which='minor', length=6, width=1)

	# Get autoscaled axes limits, and set them on secondary axis
	norm_omega_lims = np.array(ax.get_ylim())
	period_max, period_min = np.abs(Period2Omega / norm_omega_lims)

	# Get secondary axis
	ax2 = ax.twinx()
	ax2.set_ylabel(r'$\mathrm{P_{rot}}$ (d)', rotation=-90, fontsize=22, labelpad=25)
	ax2.set_yscale('log')
	ax2.set_ylim(period_max, period_min)
	ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x,pos:(f'{x:.1f}') if (x < 1) else (f'{int(x)}')))
	ax2.tick_params(labelsize=16, which='major', length=10, width=2)
	ax2.tick_params(which='minor', length=6, width=1)

	ax.set_title(plot_title, fontsize=26)
	ax.legend(loc='upper right', fontsize=14)
	plt.tight_layout()
	fig.savefig('solar_periods_log.pdf')
	plt.close()


	# Create regular plot
	fig, ax = plt.subplots(figsize=(12, 8))
	ax.scatter(solar_mags, solar_periods, s=mark_size, marker='+')
	ax.axhline(slow_period, lw=1.0, c='red', label='ZAMS Slow Rot')
	# ax.axhline(med_period, lw=1.0, c='green', label='ZAMS Med Rot')
	ax.axhline(fast_period, lw=1.0, c='blue', label='ZAMS Fast Rot')
	ax.axvline(solar_jk[0], ls='--', lw=1.0, c='#a15f43')
	ax.axvline(solar_jk[1], ls='--', lw=1.0, c='#a15f43', label=r'0.9-1.1 $M_\odot$')
	ax.set_xlabel('J - K', fontsize=22)
	ax.set_ylabel(r'$\mathrm{P_{rot}}$ (d)', fontsize=22)
	ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
	ax.tick_params(top=True, which='both', labelsize=16)
	ax.tick_params(which='major', length=10, width=2)
	ax.tick_params(which='minor', length=6, width=1)

	ax.set_title(plot_title, fontsize=26)
	ax.legend(loc='upper left', fontsize=14)
	plt.tight_layout()
	fig.savefig('solar_periods.pdf')
	plt.close()
