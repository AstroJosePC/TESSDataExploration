"""
Download target pixel files cutouts for targets in our target list
Note: It will only be download the sector 8 frames

TPF Naming notation:
    {TIC ID}_TPF_S{Sector#}C{CutoutSize}.fits
"""
from os import mkdir
from os.path import isdir

from astropy.io import ascii
from lightkurve import search_tesscut

from usefulFuncs import getOrigAps

sector = 8
cutout_size = 8

# Import target list
targets = ascii.read('DataInput/cluster_targets_tic.ecsv')
gmin = targets['G Group'].min()

# Get original apertures
group_apts = getOrigAps(targets)

# Filename template
template_fn = './TESSCuts/TPF{ticid}_S{sec}C{cutsize}.fits'

# Create folder to save TPFs
if not isdir('TESSCuts'):
    mkdir('TESSCuts')

# Download each TESS cuts for each target
for ticid, Ggroup in targets[['TIC ID', 'G Group']]:
    # Get aperture and convert to Keppler/Tess format
    apt = group_apts[Ggroup - gmin] + 1
    apt[apt == 2] += 1

    print(f'Downloading TIC ID={ticid}')
    tpf = search_tesscut(ticid, sector=sector).download(cutout_size=cutout_size)

    # Write TIC ID and pipeline aperture to HDUs
    tpf.hdu[0].header['TIC ID'] = ticid
    tpf.hdu[0].header['TICID'] = ticid
    tpf.hdu[2].data = apt

    tpf_fn = template_fn.format(ticid=ticid, sec=sector, cutsize=cutout_size)

    print(f'Saving to {tpf_fn}')
    tpf.to_fits(output_fn=tpf_fn, overwrite=True)
