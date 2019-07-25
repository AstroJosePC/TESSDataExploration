from astropy.coordinates import SkyCoord
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Column
from astroquery.mast import Tesscut
from lightkurve import TessTargetPixelFile, search_targetpixelfile, search_tesscut, TessLightCurve, TessLightCurveFile
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from random import randint

# Some constants
sector = 8
cutoutsize = 8
lc_filename = '../LightCurves/TESS_LC_{0}_SEC{1}.fits'

# List of SampleTPF files that have the aperture masks
sample_tpfs = glob('./SampleTPFs/*.fits')

# Make parallel list of target names; convert ID's into integers
sample_targets = [int(filename.split('_')[1]) for filename in sample_tpfs]

# Import target list
targets = ascii.read("../Tables/cluster_targets_tic.csv")

# This is the mags we loop over
mags = np.arange(10, 20)

# Loop over magnitude groups; e.g. 10 mag, 12 mag, 13 mag, etc.
for mag in mags:
    print('Currently producing light curves for mag=', mag)
    # All targets of current mag
    idx = (mag <= targets['G']) & (targets['G'] < mag+1)
    group_targets = targets[idx]['TIC ID']
    
    # Get sample tpfs for current group of targets
    # group_samples = [target for target in sample_targets if target in group_targets]
    group_samples = [sample_tpfs[i] for i in range(len(sample_targets)) if sample_targets[i] in group_targets]
    
    # If there are no group samples, continue with the next one
    if len(group_samples) == 0:
        print('No sample tpfs found for current group')
        continue
    
    # Assuming we got at least one sample tpf, extract their aperture masks
    apertures = []
    for sample in group_samples:
        print('Gathering aperture masks from sample TPF', sample)
        lcf = TessLightCurveFile(sample)
        current_aperture = lcf.hdu[2].data
        apertures.append(current_aperture)
    
    # Generate boolean master aperture mask 
    # I want to keep the mask pixels that are shared between at least two images. 
    # That means that if I sum the aperture masks, I will keep the mask pixels that are above 6.
    boolean_apt = np.sum(apertures, axis=0) > 6.
    
    # Now use the boolean aperture mask to export light curves of current group
    for target in group_targets:
        print('Exporting light curve file for TIC =', target)
        # Download cutouts to create light curve
        tpf = search_tesscut(target, sector=sector).download(cutout_size=cutoutsize)
        # Generate light curve using the aperture mask `new_apt`
        lc = tpf.to_lightcurve(aperture_mask=boolean_apt)
        # Export light curve into FITS file
        lc.to_fits(lc_filename.format(target, sector), overwrite=True)
        # scatter plot light curve
        # lc.scatter()
        
        