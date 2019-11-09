import numpy as np
import matplotlib.pyplot as plt
from lightkurve import MPLSTYLE
from lightkurve import open as open_tpf
from os.path import isfile, join
from os import listdir
import warnings
import re


def plot_frame(tpf, aperture=None, ax=None, savefn=None, frame=200, show_colorbar=True, **kwargs):
    """
    Function to plot TPF frame
    (optional) If aperture is 'threshold' string, a mask will be generated from
    the TPF native create_threshold_mask() function.

    :param tpf: target pixel file
    :param aperture: aperture mask if any to plot over image
    :param ax: if given, the frame will plot on this axis
    :param savefn: safe filename if want to save file
    :param frame: frame # of TPF which to plot
    :param show_colorbar: whether to show color bar
    :param kwargs: keyword arguments to pass to tpf.plot function
    """
    if not ax:
        ax = plt.subplot(projection=tpf.wcs)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        if aperture == 'threshold':
            aperture = tpf.create_threshold_mask()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        tpf.plot(ax=ax, interpolation="nearest", aperture_mask=aperture, cmap="hot",
                 scale='sqrt', frame=frame, show_colorbar=show_colorbar, **kwargs)
    with plt.style.context(MPLSTYLE):
        ax.set_ylabel('Declination')
        ax.set_xlabel('Right Ascension')
    # IF want to save
    if savefn:
        plt.gcf().savefig(savefn)


def find_tpf(ticid, sector=8, cutsize=8):
    fname = f'./TESSCuts8/TPF{ticid}_S{sector}C{cutsize}.fits'
    if isfile(fname):
        return open_tpf(fname)
    else:
        print(f'Target Pixel File w/ TICID={ticid} not found')


def find_tpfs(ticids=None, sectors=None, cutsizes=None):
    """
    Retrieve Target Pixel Files from saved cutouts
    :param ticids: list of TIC IDs of targets, string or list
    :param sectors: sectors for which to retrieve, int or list
    :param cutsizes: cut out sizes for which to retrieve, int or list
    :return: list of target pixel file objects
    """
    return [*ifind_tpfs(ticids=ticids, sectors=sectors, cutsizes=cutsizes)]


def ifind_tpfs(ticids=None, sectors=None, cutsizes=None):
    """
    Retrieve Target Pixel Files from saved cutouts; as a generator
    :param ticids: list of TIC IDs of targets, string or list
    :param sectors: sectors for which to retrieve, int or list
    :param cutsizes: cut out sizes for which to retrieve, int or list
    :return: generator of target pixel file objects
    """
    template_fn = r'^TPF{ticid}_S{sec}C{cut}.fits$'
    cuts = secs = r'\d{1,2}'

    if sectors:
        sectorsarr = np.atleast_1d(sectors)
        secs = fr'({"|".join(sectorsarr)})'
    if cutsizes:
        cutsarr = np.atleast_1d(cutsizes)
        cuts = fr'({"|".join(cutsarr)})'

    all_tpfs = listdir('TESSCuts')

    if not ticids:
        ticidsre = r'\d{8,9}'
    else:
        ticidsarr = np.atleast_1d(ticids)
        ticidsre = fr'({"|".join(ticidsarr)})'

    regexp = re.compile(template_fn.format(ticid=ticidsre, sec=secs, cut=cuts))
    for tpf in all_tpfs:
        if regexp.match(tpf):
            tpf_fn = join('.', 'TESSCuts', tpf)
            yield open_tpf(tpf_fn)
