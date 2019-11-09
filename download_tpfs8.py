"""
Download target pixel files cutouts for targets in our target list
Note: It will only be download the sector 8 frames

TPF Naming notation:
    {TIC ID}_TPF_S{Sector#}C{CutoutSize}.fits
"""
from lightkurve import search_tesscut
from astropy.io import ascii

# Parameters
sector = 8
cutout_size = 8

# Import target list
targets = ascii.read('DataInput/cluster_targets_tic.ecsv')

# Filename template
template_fn = './TESSCuts/TPF{ticid}_S{sec}C{cutsize}.fits'

# Download each TESS cuts for each target
for target in targets:
    ticid = target['TIC ID']
    
    print(f'Downloading TIC ID={ticid}')
    tpf = search_tesscut(ticid, sector=sector).download(cutout_size=cutout_size)
    tpf.hdu[0].header['TICID'] = ticid
    tpf_fn = template_fn.format(ticid=ticid, sec=sector, cutsize=cutout_size)
    
    print(f'Saving to {tpf_fn}')
    tpf.to_fits(output_fn=tpf_fn, overwrite=True)