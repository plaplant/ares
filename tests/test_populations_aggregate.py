"""

test_populations_bh.py

Author: Jordan Mirocha
Affiliation: McGill
Created on: Fri  3 Apr 2020 12:54:41 EDT

Description:

"""

import ares
import numpy as np
from ares.physics.Constants import rhodot_cgs

def test():
    pop = ares.populations.GalaxyPopulation(pop_fstar=0.1, pop_sfr_model='fcoll',
        pop_sfrd=None)

    zarr = np.arange(5, 40)

    sfrd = pop.get_sfrd(zarr) * rhodot_cgs

    zeta = pop.get_zeta_ion(6)

    assert zeta == 40, "Default parameters should yield zeta=40."

    # Supply same SFRD as function, make sure we get the same answer.
    pop3 = ares.populations.GalaxyPopulation(pop_sfr_model='sfrd-func',
        pop_sfrd=lambda z: np.interp(z, zarr, sfrd))

    assert abs(sfrd[5] - pop3.get_sfrd(zarr[5])) < 1e-8, \
        "{:.3e} {:.3e}".format(sfrd[5], pop3.SFRD(zarr[5]))

    # Make sure we get zero outside (pop_zdead, pop_zform)
    assert pop.get_sfrd(100) == 0

if __name__ == '__main__':
    test()
