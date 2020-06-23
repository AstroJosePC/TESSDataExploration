"""
Download target pixel files cutouts for targets in our target list.
The following info mentions handpicked apertures. By this I mean the grouping that I made
to assign a handpicked aperture to multiple targets of common magnitude range to reduce 
the manual aperture picking. Unless you are familiar with this you should probably 
not use this option.

There are a few paramters to take care of:

:param targets_file: Filepath to target list, assumes astropy compatible table format
:param sector: sector number to download data from; default Sector 8.
:param cutout_size: frame cutout size to get; default 8x8
:param ticid_col: Column in input table that has the targets TIC ID
:param load_handpicked: Whether to load handpicked data. This will probably be false.
:param handpicked_key: Column string that contains the grouping magnitude to assign apts
TPF Naming notation:
    {TIC ID}_TPF_S{Sector#}C{CutoutSize}.fits
"""
from os import mkdir
from os.path import isdir

from astropy.io import ascii
from lightkurve import search_tesscut

from usefulFuncs import getOrigAps

# Parameters for program
targets_file = 'DataInput/cluster_targets_tic.ecsv'
sector = 8
cutout_size = 8
ticid_col = 'TIC ID'
load_handpicked = False
handpicked_key = 'G Group'

# Filename template
template_fn = './TESSCuts/TPF{ticid}_S{sec}C{cutsize}.fits'

# Create folder to save TPFs
if not isdir('TESSCuts'):
    mkdir('TESSCuts')

# Import target list
targets = ascii.read(targets_file)

if load_handpicked:
    # Get hand picked group key
    gmin = targets[handpicked_key].min()

    # Get original apertures
    group_apts = getOrigAps(targets)

    # Download each TESS cuts for each target
    for ticid, Ggroup in targets[[ticid_col, handpicked_key]]:
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
else:
    # Download each TESS cuts for each target
    for ticid in targets[ticid_col]:
        print(f'Downloading TIC ID={ticid}')
        tpf = search_tesscut(ticid, sector=sector).download(cutout_size=cutout_size)

        # Write TIC ID and pipeline aperture to HDUs
        tpf.hdu[0].header['TIC ID'] = ticid
        tpf.hdu[0].header['TICID'] = ticid

        tpf_fn = template_fn.format(ticid=ticid, sec=sector, cutsize=cutout_size)

        print(f'Saving to {tpf_fn}')
        tpf.to_fits(output_fn=tpf_fn, overwrite=True)