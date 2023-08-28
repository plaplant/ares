"""
Behroozi et al. (2019) compilation of observed data.
"""

import os
import numpy as np
from ares.data import ARES

info = \
{
 'reference': 'Behroozi, Wechsler, Hearin, & Conroy, 2019, MNRAS, 488, 3143',
}

_input = ARES + '/umachine-data/umachine-dr1'

def get_data(field, flag=None, sources=None):
    """
    field options are 'smf', 'uvlf', 'ssfr', 'csfr', 'qf'

    .. note :: For stellar mass functions, things are trickier if we want to
        sub-divide based on star-forming vs. quiescent galaxies. This is the
        the sole purpose of the `flag` keyword argument so far. Set `flag=q`
        for quiescent or 'sf' for star-forming when field='smf' to subdivide.
    """

    data = {}
    for fn in os.listdir(_input + '/observational_constraints'):
        if not fn.endswith(field):
            continue

        if field == 'csfr':
            src = fn[0:fn.find('.')]
        else:
            src = fn[0:fn.find('_')]

        if sources is not None:
            if src not in sources:
                continue

        if src not in data:
            data[src] = {}

        # First, read-in header only to figure out what we're dealing with.
        f = open(f"{_input}/observational_constraints/{fn}", 'r')

        hdr = {}
        line = f.readline()
        while line.startswith('#'):
            name, colon, val = line[1:].partition(":")
            hdr[name] = val.strip()
            line = f.readline()
        f.close()

        # Special case: cosmic SFRD
        if hdr['type'].startswith('cosmic sfr'):
            zlo, zhi, logSFRD, errlo, errhi = \
                np.loadtxt(f"{_input}/observational_constraints/{fn}", unpack=True)

            if type(zlo) in [int, float, np.float64]:
                zarr = np.array([[zlo], [zhi]]).T
                err = np.array([[errlo], [errhi]]).T
            else:
                zarr = np.array([zlo, zhi]).T
                err = np.atleast_2d(np.array([errlo, errhi])).T

            data[src] = zarr, logSFRD, err, hdr

            continue

        if 'zlow' not in hdr:
            print(f"Dunno what to do with {src}")
            continue

        ##
        # Everything else
        zlo, zhi = float(hdr['zlow']), float(hdr['zhigh'])
        zmean = (zlo + zhi) / 2.

        Mlo, Mhi, phi, errlo, errhi = \
            np.loadtxt(f"{_input}/observational_constraints/{fn}",
            unpack=True)

        if np.any(phi < 0):
            phi_is_log = True
        else:
            phi_is_log = False

        xerr = (Mhi - Mlo) / 2.
        M = (Mhi + Mlo) / 2.

        if phi_is_log:
            yerr = np.array([errlo, errhi]).T
        else:
            yp = np.log10(phi + errhi) - np.log10(phi)
            ym = np.log10(phi) - np.log10(phi - errlo)
            yerr = np.array([ym, yp])
            if yerr.ndim == 1:
                yerr = np.atleast_2d(yerr).T

        x = np.array([Mlo, Mhi]).T

        ##
        # Modifications for smf with flag='sf' or 'q'
        # Retrieve quenched fraction in this case
        if (field == 'smf') and (flag is not None):
            data_qf = get_data(field='qf', sources=src)

            # First, check that mass range is the same
            assert np.all(x == data_qf[src][(zlo, zhi)][0]), \
                f"Mismatch in {src} mass bins for SMF flag={flag}!"



        ##
        # Pack up and move on
        data[src][(zlo, zhi)] = \
            np.atleast_2d(x), phi, np.atleast_2d(yerr), hdr

    ##
    # Done
    return data

def get_results(field, method=None):
    """
    Now we're hunting for actual modeling results.
    """

    if field != 'smhm':
        raise NotImplemented('help')

    subdir = f'{_input}/data/{field}'
    if method is not None:
        subdir += f'/{method}'

    data = {}
    for fn in os.listdir(subdir):
        if not fn.startswith(f'{field}_a'):
            continue

        ##
        # Convention here is to put scale factor in filename, e.g.,
        # `field`_a<0.12345>.dat
        astr = fn[fn.rfind('_')+2:fn.rfind('.')]
        a = float(astr)
        z = (1. / a) - 1
        data[z] = {}

        # Load
        _data = np.loadtxt(f"{subdir}/{fn}", unpack=True)

        with open(f"{subdir}/{fn}", 'r') as f:
            _cols = f.readline()[1:].split()

        # Get rid of column numbers embedded in header
        cols = [col[0:col.rfind('(')] for col in _cols]

        data[z]['cols'] = cols
        data[z]['data'] = _data

    return data
