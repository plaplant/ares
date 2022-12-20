"""

SynthesisModelToy.py

Author: Jordan Mirocha
Affiliation: McGill
Created on: Mon 20 May 2019 15:08:49 EDT

Description:

"""

import numpy as np
from ..util.Stats import bin_c2e
from ..util.ReadData import read_lit
from functools import cached_property
from .SynthesisModel import SynthesisModelBase
from ..util.ParameterFile import ParameterFile
from ..physics.Constants import c, h_p, erg_per_ev, cm_per_ang, ev_per_hz, E_LL

class SynthesisModelToy(SynthesisModelBase):
    def __init__(self, **kwargs):
        SynthesisModelBase.__init__(self, **kwargs)
        #self.pf = ParameterFile(**kwargs)

        #self.Emin = self.pf['source_Emin']
        #self.Emax = self.pf['source_Emax']

    @property
    def tab_energies_c(self):
        if not hasattr(self, '_energies'):
            if (self.pf['source_wavelengths'] is not None) or \
               (self.pf['source_lmin'] is not None):
                self._energies = h_p * c / self.tab_waves_c / cm_per_ang \
                    / erg_per_ev
            else:
                dE = self.pf['source_dE']
                Emin = self.pf['source_Emin']
                Emax = self.pf['source_Emax']

                self._energies = np.arange(Emin, Emax+dE, dE)[-1::-1]

            # Should be in descending order
            assert np.all(np.diff(self._energies) < 0)

        return self._energies

    @property
    def tab_waves_c(self):
        if not hasattr(self, '_wavelengths'):
            if self.pf['source_wavelengths'] is not None:
                self._wavelengths = self.pf['source_wavelengths']
            elif self.pf['source_lmin'] is not None:
                self._wavelengths = np.arange(self.pf['source_lmin'],
                    self.pf['source_lmax']+self.pf['source_dlam'],
                    self.pf['source_dlam'])
            else:
                self._wavelengths = h_p * c / self.tab_energies_c / erg_per_ev \
                    / cm_per_ang

            if (self._wavelengths.max() < 2e3):
                raise ValueError('Wavelengths all FUV. This is generally not a good idea!')

            # Should be in ascending order
            assert np.all(np.diff(self._wavelengths) > 0)

        return self._wavelengths

    @cached_property
    def tab_waves_e(self):
        self._waves_e = bin_c2e(self.tab_waves_c)
        return self._waves_e

    @property
    def frequencies(self):
        if not hasattr(self, '_frequencies'):
            self._frequencies = c / self.tab_waves_c / cm_per_ang
        return self._frequencies

    @property
    def tab_t(self):
        if not hasattr(self, '_times'):
            if self.pf['source_times'] is not None:
                self._times = self.pf['source_times']
            else:
                # Standard SPS time gridding
                self._times = 10**np.arange(0, 4.1, 0.1)
        return self._times

    @cached_property
    def tab_dE(self):
        self._dE = np.diff(self.tab_energies_e)
        return self._dE

    @cached_property
    def tab_dndE(self):
        self._dndE = np.abs(np.diff(self.tab_freq_e) / np.diff(self.tab_energies_e))
        return self._dndE

    @cached_property
    def tab_dwdn(self):
        self._dwdn = self.tab_waves_c**2 / (c * 1e8)
        return self._dwdn

    @property
    def _Star(self):
        if not hasattr(self, '_Star_'):

            #assert self.pf['source_ssp'], "Must set source_ssp=True!"

            from ..sources import StarQS, Star
            kw = self.pf.copy()
            kw['source_sed'] = self.pf['source_toysps_method']

            if self.pf['source_toysps_method'] == 'schaerer2002':
                self._Star_ = StarQS(cosm=self.cosm, **kw)
            else:
                self._Star_ = Star(cosm=self.cosm, **kw)

        return self._Star_

    def _Spectrum(self, t, wave=1600.):

        if self.pf['source_toysps_method'] == 0:
            beta = self.pf["source_toysps_beta"]
            _norm = self.pf["source_toysps_norm"]
            gamma = self.pf["source_toysps_gamma"]
            delta = self.pf["source_toysps_delta"]
            alpha = self.pf["source_toysps_alpha"]
            trise = self.pf['source_toysps_trise']
            _t0 = self.pf['source_toysps_t0']
            lmin = self.pf['source_toysps_lmin']

            ok = wave >= lmin

            # Normalization of each wavelength is set by UV slope
            norm = _norm * (wave / 1600.)**beta

            # Assume that all wavelengths initially decline as a power-law
            # with the same index
            _gamma_ = gamma * (wave / 1600.)**delta
            pl_decay = (t / 1.)**_gamma_

            # Assume an exponential decay at some critical (wavelength-dependent)
            # timescale.
            t0 = _t0 * (wave / 1600.)**alpha
            exp_decay = np.exp(-t / t0)

            # Put it all together.
            spec = norm * pl_decay * exp_decay
            spec[ok==0] = 0

            # Assume log-linear at t < trise
            if t < trise:
                spec *= (t / trise)**1.5

        # Power-law, time-independent spectrum.
        elif self.pf['source_toysps_method'] == 1:
            beta = self.pf["source_toysps_beta"]
            _norm = self.pf["source_toysps_norm"]
            lmin = self.pf['source_toysps_lmin']

            # Normalization of each wavelength is set by UV slope
            norm = _norm * (wave / 1600.)**beta
            ok = wave >= lmin
            spec = norm
            spec[ok==0] = 0

        elif type(self.pf['source_toysps_method']) == str:

            is_on = t < (self._Star.lifetime / 1e6) \
                or (not self.pf['source_ssp'])

            if is_on:
                # This is normalized to Q in each sub-band.
                E = h_p * c / (wave * cm_per_ang) / erg_per_ev
                spec = self._Star.Spectrum(E)
                # Right here, `spec` integrates to unity over relevant bands.

                mass = self._Star.pf['source_mass']

                # Here, we need to re-normalize again to Q, then switch to
                # units of erg/s/A/Msun (instead of eV^-1)

                # _StarQS.norm_ contains normalization factors that,
                # when multiplied by a blackbody in the relevant spectral range,
                # will yield Q * <Eband> * erg_per_ev

                # Make an array like `E` with relevant _StarQS.norm_ values
                if type(E) not in [list, tuple, np.ndarray]:
                    E = np.array([E])
                    spec = np.array([spec])

                if self.pf['source_toysps_method'] == 'schaerer2002':

                    # Why isn't this all handled in StarQS?
                    # Caused by new wavelength gridding?
                    bands = self._Star.bands

                    norms = []
                    for i, nrg in enumerate(E):
                        if nrg < bands[0][1]:
                            norms.append(self._Star.norm_[0])
                        elif bands[1][0] <= nrg < bands[1][1]:
                            norms.append(self._Star.norm_[1])
                        elif bands[2][0] <= nrg < bands[2][1]:
                            norms.append(self._Star.norm_[2])
                        else:
                            norms.append(self._Star.norm_[2])

                    norms = np.array(norms)
                    spec *= norms
                elif self.pf['source_toysps_method'] == 'bb':
                    spec *= self._Star.Lbol0
                else:
                    raise NotImplemented('help')

                # Now, [spec] = erg/s/eV

                # Do last unit conversion to get erg/s/A/Msun
                spec *= ev_per_hz / self.dwdn

                if self.pf['source_ssp']:
                    raise NotImplemented('Set source_ssp=False for now.')
                    spec /= mass
                else:
                    spec /= (mass / self._Star.lifetime)

            else:
                spec = 0.

        else:
            raise NotImpemented('Do not recognize source_toysps_method={}'.format(
                self.pf['source_toysps_method']
            ))

        return spec

    @property
    def tab_sed(self):
        """
        Units of erg / s / A / Msun
        """
        if not hasattr(self, '_data'):
            self._data = np.zeros((self.tab_waves_c.size, self.tab_t.size))
            for i, t in enumerate(self.tab_t):
                self._data[:,i] = self._Spectrum(t, wave=self.tab_waves_c)

            self._add_nebular_emission()

        return self._data
