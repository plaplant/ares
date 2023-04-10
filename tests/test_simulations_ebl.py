"""

test_simulations_ebl.py

Author: Jordan Mirocha
Affiliation: JPL / Caltech
Created on: Sun Apr  9 15:30:19 PDT 2023

Description:

"""

import ares
import numpy as np

def test():
    pars = ares.util.ParameterBundle('mirocha2023:centrals_sf')
    pars.update(ares.util.ParameterBundle('testing:galaxies'))

    pop = ares.populations.GalaxyPopulation(**pars)

    P2h = pop.get_ps_2h(z=pop.halos.tab_z[0], k=1, wave1=1600., wave2=1600)
    Pshot = pop.get_ps_shot(z=pop.halos.tab_z[0], k=1, wave1=1600., wave2=1600)

    assert P2h < Pshot

    Pshot2 = pop.get_ps_shot(z=pop.halos.tab_z[0], k=10, wave1=1600., wave2=1600)
    assert Pshot == Pshot2

    # Shot noise in P(k) sense larger or smaller as z increases?
    Pshot3 = pop.get_ps_shot(z=pop.halos.tab_z[1], k=1, wave1=1600., wave2=1600)
    #assert Pshot3 > Pshot2

    # Check that fluctuations and mean are order unity in SI units


"""
Things to test:

Qualitative dependencies of shot noise and 2h power on redshift, wavelength,
sensitivity to Mmin.

"""

if __name__ == '__main__':
    test()
