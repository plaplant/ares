# Thanks, Jason Sun, for most of this!

import os
import re
import time
import pickle
from types import FunctionType, MethodType

import numpy as np
import scipy.special as sp
from scipy.integrate import quad

from ..data import ARES
from ..util.ProgressBar import ProgressBar
from .HaloMassFunction import HaloMassFunction
from .Constants import rho_cgs, c, cm_per_mpc
from ..util.Math import get_cf_from_ps_tab, get_cf_from_ps_func

try:
    import h5py
except ImportError:
    pass

try:
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
except ImportError:
    rank = 0
    size = 1

four_pi = 4 * np.pi

class HaloModel(HaloMassFunction):

    def get_concentration(self, z, Mh):
        """
        Get halo concentration from named concentration-mass-relation (CMR).

        Parameters
        ----------
        z : int, float
            Redshift
        Mh : int, float, numpy.ndarray
            Halo mass [Msun].
        return_Rs : bool
            If True, return a tuple containing (concentration, Rvir / c),
            otherwise just returns concentration.

        """
        if self.pf['halo_cmr'] == 'duffy':
            return self._get_cmr_duffy(z, Mh)
        elif self.pf['halo_cmr'] == 'zehavi':
            return self._get_cmr_zehavi(z, Mh)
        else:
            raise NotImplemented('help!')

    def _get_cmr_duffy(self, z, Mh):
        """
        Concentration-mass relation from Duffy et al. (2008).

        .. note :: This is from Table 1, row 4, which is for z=0-2.
        """
        c = 6.71 * (Mh / (2e12)) ** -0.091 * (1 + z) ** -0.44
        return c

    def _get_cmr_zehavi(self, z, Mh):
        c = ((m / 1.5e13) ** -0.13) * 9.0 / (1 + z)
        return c

    def get_Rvir_from_Mh(self, Mh):
        return (3. * Mh / (4. * np.pi * self.pf['halo_delta'] \
            * self.cosm.mean_density0)) ** (1. / 3.)

    def _dc_nfw(self, c):
        return c** 3. / (4. * np.pi) / (np.log(1 + c) - c / (1 + c))

    def get_rho_nfw(self, z, Mh, r, truncate=True):

        con = self.get_concentration(z, Mh)
        rvir = self.get_Rvir_from_Mh(Mh)
        r_s = rvir / con

        x = r / r_s
        rn = x / con

        small_sc = rn <= 1 if truncate else True

        if np.iterable(x):
            result = np.zeros_like(x)
            result[small_sc==1] = (self._dc_nfw(con) \
                / (con * r_s)**3 / (x * (1 + x)**2))[rn <= 1]

            return result
        else:
            if small_sc:
                return self._dc_nfw(c) / (con * r_s) ** 3 / (x * (1 + x) ** 2)
            else:
                return 0.0

    def get_u_nfw(self, z, Mh, k):
        """
        Normalized Fourier Transform of an NFW profile.

        ..note:: This is Equation 81 from Cooray & Sheth (2002).

        Parameters
        ----------
        z : int, float
            Redshift
        Mh : int, float, numpy.ndarray
            Halo mass [Msun].
        k : int, float
            Wavenumber

        """
        con = self.get_concentration(z, Mh)
        rvir = self.get_Rvir_from_Mh(Mh)
        r_s = rvir / con

        K = k * r_s

        asi, ac = sp.sici((1 + con) * K)
        bs, bc = sp.sici(K)

        # The extra factor of np.log(1 + c) - c / (1 + c)) comes in because
        # there's really a normalization factor of 4 pi rho_s r_s^3 / m,
        # and m = 4 pi rho_s r_s^3 * the log term
        norm = 1. / (np.log(1 + con) - con / (1 + con))

        return norm * (np.sin(K) * (asi - bs) - np.sin(con * K) / ((1 + con) * K) \
            + np.cos(K) * (ac - bc))

    @property
    def tab_u_nfw(self):
        if not hasattr(self, '_tab_u_nfw'):
            fn = os.path.join(
                ARES, "input", "halos", self.tab_prefix_prof() + ".hdf5"
            )

            if os.path.exists(fn):
                with h5py.File(fn, 'r') as f:
                    self._tab_u_nfw = np.array(f[('tab_u_nfw')])

                print(f"# Loaded {fn}.")
            else:
                self._tab_u_nfw = None
                print(f"# Did not find {fn}. Will generate u_nfw on the fly.")

        return self._tab_u_nfw

    @property
    def tab_Sigma_nfw(self):
        if not hasattr(self, '_tab_Sigma_nfw'):

            fn = f"{ARES}/input/halos/{self.tab_prefix_surf()}.hdf5"

            if os.path.exists(fn):
                with h5py.File(fn, 'r') as f:
                    self._tab_Sigma_nfw = np.array(f[('tab_Sigma_nfw')])

                print(f"# Loaded {fn}.")
            else:
                self._tab_Sigma_nfw = None
                print(f"# Did not find {fn}.")

        return self._tab_Sigma_nfw

    def get_u_isl(self, z, Mh, k, rmax=1e2):
        """
        Normalized Fourier transform of an r^-2 profile.

        Parameters
        ----------
        z : int, float
            Redshift
        Mh : int, float, numpy.ndarray
            Halo mass [Msun].
        k : int, float
            Wavenumber
        rmax : int, float
            Effective horizon. Distance a photon can travel between
            Ly-beta and Ly-alpha.

        """

        asi, aco = sp.sici(rmax * k)

        return asi / rmax / k

    def get_u_isl_exp(self, z, Mh, k, rmax=1e2, rstar=10):
        return np.arctan(rstar * k) / rstar / k

    def get_u_exp(self, z, Mh, k, rmax=1e2):
        rs = 1.

        L0 = (Mh / 1e11)**1.
        c = rmax / rs

        kappa = k * rs

        norm = rmax / rs**3

        return norm / (1. + kappa**2)**2.

    def get_u_cgm_rahmati(self, z, Mh, k):
        rstar = 0.0025
        return np.arctan((rstar * k) ** 0.75) / (rstar * k) ** 0.75

    def get_u_cgm_steidel(self, z, Mh, k):
        rstar = 0.2
        return np.arctan((rstar * k) ** 0.85) / (rstar * k) ** 0.85

    def FluxProfile(self, r, m, z, lc=False):
        raise NotImplemented('need to fix this')
        return m * self.ModulationFactor(z, r=r, lc=lc) / (4. * np.pi * r**2)

    #@RadialProfile.setter
    #def RadialProfile(self, value):
    #    pass

    def FluxProfileFT(self, k, m, z, lc=False):
        raise NotImplemented('need to fix this')
        _numerator = lambda r: 4. * np.pi * r**2 * np.sin(k * r) / (k * r) \
            * self.FluxProfile(r, m, z, lc=lc)
        _denominator = lambda r: 4. * np.pi * r**2 *\
            self.FluxProfile(r, m, z, lc=lc)
        _r_LW = 97.39 * self.ScalingFactor(z)
        temp = quad(_numerator, 0., _r_LW)[0] / quad(_denominator, 0., _r_LW)[0]
        return temp

    def ScalingFactor(self, z):
        return (self.cosm.h70 / 0.7)**-1 \
            * (self.cosm.omega_m_0 / 0.27)**-0.5 * ((1. + z) / 21.)**-0.5

    def ModulationFactor(self, z0, z=None, r=None, lc=False):
        """
        Return the modulation factor as a function of redshift or comoving distance
        - Reference: Ahn et al. 2009
        :param z0: source redshift
        :param z: the redshift (whose LW intensity is) of interest
        :param r: the distance from the source in cMpc
        :lc: True or False, including the light cone effect
        :return:
        """
        if z != None and r == None:
            r_comov = self.cosm.ComovingRadialDistance(z0, z)
        elif z == None and r != None:
            r_comov = r
        else:
            raise ValueError('Must specify either "z" or "r".')
        alpha = self.ScalingFactor(z0)
        _a = 0.167
        r_star = c * _a * self.cosm.HubbleTime(z0) * (1.+z0) / cm_per_mpc
        ans = np.maximum(1.7 * np.exp(-(r_comov / 116.29 / alpha)**0.68) - 0.7, 0.0)
        if lc == True:
            ans *= np.exp(-r/r_star)
        return ans

    def _get_ps_integrals(self, k, iz, prof1, prof2, lum1, lum2, mmin1, mmin2,
        focc1, focc2, term):
        """
        Compute integrals over profile, weighted by bias, dndm, etc.,
        needed for halo model.

        .. note :: This is really just a wrapper around _integrate_over_prof,
            that handles the fact that `k` can be a number or an array.

        """

        if type(k) == np.ndarray:
            integ1 = []; integ2 = []
            for _k in k:
                _integ1, _integ2 = self._integrate_over_prof(_k, iz,
                    prof1, prof2, lum1, lum2, mmin1, mmin2, focc1, focc2, term)
                integ1.append(_integ1)
                integ2.append(_integ2)

            integ1 = np.array(integ1)
            integ2 = np.array(integ2)
        else:
            integ1, integ2 = self._integrate_over_prof(k, iz,
                prof1, prof2, lum1, lum2, mmin1, mmin2, focc1, focc2, term)

        return integ1, integ2

    def _integrate_over_prof(self, k, iz, prof1, prof2, lum1, lum2, mmin1,
        mmin2, focc1, focc2, term):
        """
        Compute integrals over profile, weighted by bias, dndm, etc.,
        needed for halo model.
        """

        if type(prof1) in [FunctionType, MethodType]:
            p1 = np.abs([prof1(self.tab_z[iz], M, k) for M in self.tab_M])
        else:
            p1 = np.abs([np.interp(k, self.tab_k, prof1[iz,iM,:]) \
                for iM in np.arange(self.tab_M.size)])

        if type(prof2) in [FunctionType, MethodType]:
            p2 = np.abs([prof2(self.tab_z[iz], M, k) for M in self.tab_M])
        else:
            p2 = np.abs([np.interp(k, self.tab_k, prof2[iz,iM,:]) \
                for iM in np.arange(self.tab_M.size)])

        bias = self.tab_bias[iz]
        rho_bar = self.cosm.rho_m_z0 * rho_cgs
        dndlnm = self.tab_dndlnm[iz]

        if (mmin1 is None) and (lum1 is None):
            fcoll1 = 1.

            # Small halo correction. Make use of Cooray & Sheth Eq. 71
            _integrand = dndlnm * (self.tab_M / rho_bar) * bias
            corr1 = 1. - np.trapz(_integrand, x=np.log(self.tab_M))
        elif lum1 is not None:
            corr1 = 0.0
            fcoll1 = 1.
        else:
            fcoll1 = self.tab_fcoll[iz,np.argmin(np.abs(mmin1-self.tab_M))]
            corr1 = 0.0

        if (mmin2 is None) and (lum2 is None):
            fcoll2 = 1.#self.mgtm[iz,0] / rho_bar
            _integrand = dndlnm * (self.tab_M / rho_bar) * bias
            corr2 = 1. - np.trapz(_integrand, x=np.log(self.tab_M))
        elif lum2 is not None:
            corr2 = 0.0
            fcoll2 = 1.
        else:
            fcoll2 = self.fcoll_2d(z, np.log10(Mmin_2))#self.fcoll_Tmin[iz]
            corr2 = 0.0

        ok = self.tab_fcoll[iz] > 0

        # If luminosities passed, then we must cancel out a factor of halo
        # mass that generally normalizes the integrand.
        if lum1 is None:
            weight1 = self.tab_M
            norm1 = rho_bar * fcoll1
        else:
            weight1 = lum1
            norm1 = 1.

        if lum2 is None:
            weight2 = self.tab_M
            norm2 = rho_bar * fcoll2
        else:
            weight2 = lum2
            norm2 = 1.

        ##
        # Are we doing the 1-h or 2-h term?
        if term == 1:
            integrand = dndlnm * focc1 * weight1 * weight2 \
                * p1 * p2 / norm1 / norm2

            result = np.trapz(integrand[ok==1], x=np.log(self.tab_M[ok==1]))

            return result, None

        elif term == 2:
            integrand1 = dndlnm * focc1 * weight1 * p1 * bias / norm1
            integrand2 = dndlnm * focc2 * weight2 * p2 * bias / norm2

            integral1 = np.trapz(integrand1[ok==1], x=np.log(self.tab_M[ok==1]),
                axis=0)
            integral2 = np.trapz(integrand2[ok==1], x=np.log(self.tab_M[ok==1]),
                axis=0)

            return integral1 + corr1, integral2 + corr2

        else:
            raise NotImplemented('dunno man')

    def _prep_for_ps(self, z, k, prof1, prof2, ztol):
        """
        Basic prep: fill prof1=None or prof2=None with defaults, determine
        the index of the requested redshift in our lookup tables.
        """

        iz = np.argmin(np.abs(z - self.tab_z))

        if abs(self.tab_z[iz] - z) > ztol:
            raise ValueError('Requested z={} not in grid (ztol={}).'.format(z,
                ztol))

        if prof1 is None:
            prof1 = self.get_u_nfw
        if prof2 is None:
            prof2 = prof1

        if k is None:
            k = self.tab_k_lin

        return iz, k, prof1, prof2

    def _get_ps_lin(self, k, iz):
        """
        Return linear matter power spectrum for requested wavenumber `k`.

        .. note :: Assumes we already know the index of the redshift of interest
            in our lookup tables, `iz`.

        """
        if k is None:
            k = self.tab_k_lin
            ps_lin = self.tab_ps_lin[iz]
        else:
            ps_lin = np.exp(np.interp(np.log(k), np.log(self.tab_k_lin),
                np.log(self.tab_ps_lin[iz])))

        return ps_lin

    def get_ps_1h(self, z, k=None, prof1=None, prof2=None, lum1=None, lum2=None,
        mmin1=None, mmin2=None, focc1=1, focc2=1, ztol=1e-3):
        """
        Compute 1-halo term of power spectrum.
        """

        iz, k, prof1, prof2 = self._prep_for_ps(z, k, prof1, prof2, ztol)

        integ1, none = self._get_ps_integrals(k, iz, prof1, prof2,
            lum1, lum2, mmin1, mmin2, focc1, focc2, term=1)

        return integ1

    def get_ps_2h(self, z, k=None, prof1=None, prof2=None, lum1=None, lum2=None,
        mmin1=None, mmin2=None, focc1=1, focc2=1, ztol=1e-3):
        """
        Get 2-halo term of power spectrum.
        """

        iz, k, prof1, prof2 = self._prep_for_ps(z, k, prof1, prof2, ztol)

        ps_lin = self._get_ps_lin(k, iz)

        # Cannot return unmodified P_lin unless no L's or Mmin's passed!
        if self.pf['halo_ps_linear']:
            if (lum1 is None) and (lum2 is None) and (mmin1 is None) and (mmin2 is None):
                return ps_lin

        integ1, integ2 = self._get_ps_integrals(k, iz, prof1, prof2,
            lum1, lum2, mmin1, mmin2, focc1, focc2, term=2)

        ps = integ1 * integ2 * ps_lin

        return ps

    def get_ps_shot(self, z, k=None, lum1=None, lum2=None, mmin1=None, mmin2=None,
        focc1=1, focc2=1, ztol=1e-3):
        """
        Compute the shot noise term quickly.
        """

        iz, k, _prof1_, _prof2_ = self._prep_for_ps(z, k, None, None, ztol)

        dndlnm = self.tab_dndlnm[iz]
        integrand = dndlnm * focc1 * lum1 * lum2
        shot = np.trapz(integrand, x=np.log(self.tab_M), axis=0)

        return shot

    def get_ps_mm(self, z, k=None, prof1=None, prof2=None, lum1=None, lum2=None,
        mmin1=None, mmin2=None, ztol=1e-3):
        """
        Return total power spectrum as sum of 1h and 2h terms.
        """

        if self.pf['halo_ps_linear']:
            ps_1h = 0
        else:
            ps_1h = self.get_ps_1h(z, k, prof1, prof2, lum1, lum2, mmin1, mmin2, ztol)

        ps_2h = self.get_ps_2h(z, k, prof1, prof2, lum1, lum2, mmin1, mmin2, ztol)

        return ps_1h + ps_2h

    def get_cf_mm(self, z, R=None, load=True, ztol=1e-2):
        """
        Compute the correlation function of the matter power spectrum.

        Parameters
        ----------
        z : int, float
            Redshift of interest.
        R : int, float, np.ndarray
            Scale(s) of interest. If not supplied, will default to self.tab_R.

        Returns
        -------
        Tuple containing (R, CF). If `R` is not supplied by user, will default
        to 1 / self.tab_k.

        """

        ##
        # Load from table if one exists.
        ##
        if self.pf['halo_ps_load'] and load and (not self.pf['halo_ps_linear']):
            iz = np.argmin(np.abs(z - self.tab_z))

            assert abs(z - self.tab_z[iz]) < ztol, \
                'Supplied redshift (%g) not in table!' % z

            if R is not None:
                if len(R) == len(self.tab_R):
                    if np.allclose(R, self.tab_R):
                        return self.tab_cf_mm[iz]
                    else:
                        return np.interp(R, self.tab_R, self.tab_cf_mm[iz])

        ##
        # Otherwise, compute PS then inverse transform to obtain CF.
        ##
        if self.pf['use_mcfit']:
            k = self.tab_k_lin if self.pf['halo_ps_linear'] else self.tab_k
            Pofk = self.get_ps_mm(z, k)
            _R_, _cf_ = get_cf_from_ps_tab(k, Pofk)

            if R is not None:
                cf = np.interp(R, _R_, _cf_)
            else:
                R = _R_
                cf = _cf_

        else:
            if R is None:
                R = self.tab_R

            cf = get_cf_from_ps_func(R, lambda kk: self.get_ps_mm(z, kk))

        return R, cf

    @property
    def tab_k(self):
        """
        k-vector constructed from hps parameters.
        """
        if not hasattr(self, '_tab_k'):
            dlogk = self.pf['halo_dlnk']
            kmi, kma = self.pf['halo_lnk_min'], self.pf['halo_lnk_max']
            logk = np.arange(kmi, kma+dlogk, dlogk)
            self._tab_k = np.exp(logk)

        return self._tab_k

    @tab_k.setter
    def tab_k(self, value):
        self._tab_k = value

    @property
    def tab_R(self):
        """
        R-vector constructed from mps parameters.
        """
        if not hasattr(self, '_tab_R'):
            dlogR = self.pf['halo_dlnR']
            Rmi, Rma = self.pf['halo_lnR_min'], self.pf['halo_lnR_max']
            logR = np.arange(Rmi, Rma+dlogR, dlogR)
            self._tab_R = np.exp(logR)

        return self._tab_R

    #@property
    #def tab_z_ps(self):
    #    """
    #    Redshift array -- different than HMF redshifts!
    #    """
    #    if not hasattr(self, '_tab_z_ps'):
    #        zmin = self.pf['halo_zmin']
    #        zmax = self.pf['halo_zmax']
    #        dz = self.pf['halo_dz']

    #        Nz = int(round(((zmax - zmin) / dz) + 1, 1))
    #        self._tab_z_ps = np.linspace(zmin, zmax, Nz)

    #    return self._tab_z_ps

    #@tab_z_ps.setter
    #def tab_z_ps(self, value):
    #    self._tab_z_ps = value

    @tab_R.setter
    def tab_R(self, value):
        self._tab_R = value

        print('Setting R attribute. Should verify it matches PS.')

    def __getattr__(self, name):

        if hasattr(HaloMassFunction, name):
            return HaloMassFunction.__dict__[name].__get__(self, HaloMassFunction)

        if (name[0] == '_'):
            raise AttributeError('This will get caught. Don\'t worry! {}'.format(name))

        if name not in self.__dict__.keys():
            if self.pf['halo_mf_load']:
                self._load_hmf()
            else:
                # Can generate on the fly!
                if name == 'tab_MAR':
                    self.generate_mar()
                else:
                    self.generate_hmf(save_MAR=False)

            if name not in self.__dict__.keys():
                self._load_ps()

        return self.__dict__[name]

    def _load_ps(self, suffix='hdf5'):
        """ Load table from HDF5 or binary. """

        if self.pf['halo_ps_linear']:
            self._tab_ps_mm = np.zeros((self.tab_z.size, self.tab_k.size))
            self._tab_cf_mm = np.zeros((self.tab_z.size, self.tab_R.size))
            for i, _z_ in enumerate(self.tab_z):
                iz = np.argmin(np.abs(_z_ - self.tab_z))
                self._tab_ps_mm[i,:] = self._get_ps_lin(self.tab_k, iz)
                R, cf = get_cf_from_ps_tab(self.tab_k, self._tab_ps_mm[i,:])
                self._tab_cf_mm[i,:] = np.interp(np.log(R), np.log(self.tab_R),
                    cf)

            return

        fn = f"{self.tab_prefix_ps()}.{suffix}"
        fn = os.path.join(ARES, "input", "halos", fn)

        if re.search('.hdf5', fn) or re.search('.h5', fn):
            f = h5py.File(fn, 'r')
            self.tab_z = f['tab_z_ps'].value
            self.tab_R = f['tab_R'].value
            self.tab_k = f['tab_k'].value
            self.tab_ps_mm = f['tab_ps_mm'].value
            self.tab_cf_mm = f['tab_cf_mm'].value
            f.close()
        elif re.search('.pkl', fn):
            f = open(fn, 'rb')
            self.tab_z = pickle.load(f)
            self.tab_R = pickle.load(f)
            self.tab_k = pickle.load(f)
            self.tab_ps_mm = pickle.load(f)
            self.tab_cf_mm = pickle.load(f)
            f.close()
        else:
            raise IOError('Unrecognized format for halo_table.')

    def tab_prefix_prof(self):
        M1, M2 = self.pf['halo_logMmin'], self.pf['halo_logMmax']
        z1, z2 = self.pf['halo_zmin'], self.pf['halo_zmax']

        dlogk = self.pf['halo_dlnk']
        kmi, kma = self.pf['halo_lnk_min'], self.pf['halo_lnk_max']

        logMsize = (self.pf['halo_logMmax'] - self.pf['halo_logMmin']) \
            / self.pf['halo_dlogM']
        zsize = ((self.pf['halo_zmax'] - self.pf['halo_zmin']) \
            / self.pf['halo_dz']) + 1

        assert logMsize % 1 == 0
        logMsize = int(logMsize)
        assert zsize % 1 == 0
        zsize = int(round(zsize, 1))

        # Should probably save NFW information etc. too
        return 'halo_prof_%s_%s_logM_%s_%i-%i_z_%s_%i-%i_lnk_%.1f-%.1f_dlnk_%.3f' \
            % (self.pf['halo_profile'], self.pf['halo_cmr'],
                logMsize, M1, M2, zsize, z1, z2, kmi, kma, dlogk)

    def tab_prefix_surf(self):
        M1, M2 = self.pf['halo_logMmin'], self.pf['halo_logMmax']

        zstr = self.get_table_zstr()

        # Hard-coded for now, change this.
        Rmi, Rma = -3, 1
        dlogR = 0.05

        logMsize = (self.pf['halo_logMmax'] - self.pf['halo_logMmin']) \
            / self.pf['halo_dlogM']

        Rsize = (Rma - Rmi) / dlogR

        assert logMsize % 1 == 0
        logMsize = int(logMsize)
        assert Rsize % 1 == 0
        Rsize = int(round(Rsize, 1))

        return 'halo_surf_%s_logM_%i_%i-%i_%s_logR_%.1f-%.1f_dlnR_%.3f' \
            % (self.pf['halo_cmr'],
                logMsize, M1, M2, zstr, Rmi, Rma, dlogR)



    def tab_prefix_ps(self, with_size=True):
        """
        What should we name this table?

        Convention:
        ps_FIT_logM_nM_logMmin_logMmax_z_nz_

        Read:
        halo mass function using FIT form of the mass function
        using nM mass points between logMmin and logMmax
        using nz redshift points between zmin and zmax

        """

        M1, M2 = self.pf['halo_logMmin'], self.pf['halo_logMmax']
        z1, z2 = self.pf['halo_zmin'], self.pf['halo_zmax']

        dlogk = self.pf['halo_dlnk']
        kmi, kma = self.pf['halo_lnk_min'], self.pf['halo_lnk_max']
        #logk = np.arange(kmi, kma+dlogk, dlogk)
        #karr = np.exp(logk)

        dlogR = self.pf['halo_dlnR']
        Rmi, Rma = self.pf['halo_lnR_min'], self.pf['halo_lnR_max']
        #logR = np.arange(np.log(Rmi), np.log(Rma)+dlogR, dlogR)
        #Rarr = np.exp(logR)


        if with_size:
            logMsize = (self.pf['halo_logMmax'] - self.pf['halo_logMmin']) \
                / self.pf['halo_dlogM']
            zsize = ((self.pf['halo_zmax'] - self.pf['halo_zmin']) \
                / self.pf['halo_dz']) + 1

            assert logMsize % 1 == 0
            logMsize = int(logMsize)
            assert zsize % 1 == 0
            zsize = int(round(zsize, 1))

            # Should probably save NFW information etc. too
            return 'halo_ps_%s_logM_%s_%i-%i_z_%s_%i-%i_lnR_%.1f-%.1f_dlnR_%.3f_lnk_%.1f-%.1f_dlnk_%.3f' \
                % (self.pf['halo_mf'], logMsize, M1, M2, zsize, z1, z2,
                   Rmi, Rma, dlogR, kmi, kma, dlogk)
        else:
            raise NotImplementedError('help')

    def tab_prefix_ps_check(self, with_size=True):
        """
        A version of the prefix to be used only for checkpointing.

        This just means take the full prefix and hack out the bit with the
        redshift interval.
        """

        prefix = self.tab_prefix_ps(with_size)

        iz = prefix.find('_z_')
        iR = prefix.find('_lnR_')

        return prefix[0:iz] + prefix[iR:]

    @property
    def tab_ps_mm(self):
        if not hasattr(self, '_tab_ps_mm'):
            self._load_ps()
        return self._tab_ps_mm

    @tab_ps_mm.setter
    def tab_ps_mm(self, value):
        self._tab_ps_mm = value

    @property
    def tab_cf_mm(self):
        if not hasattr(self, '_tab_cf_mm'):
            ps = self.tab_ps_mm
        return self._tab_cf_mm

    @tab_cf_mm.setter
    def tab_cf_mm(self, value):
        self._tab_cf_mm = value

    def TabulatePS(self, clobber=False, checkpoint=True, **ftkwargs):
        """
        Tabulate the matter power spectrum as a function of redshift and k.
        """

        pb = ProgressBar(len(self.tab_z), 'ps_dd(z)')
        pb.start()

        # Lists to store any checkpoints that are found
        _z = []
        _ps = []
        _cf = []
        if checkpoint:
            if (not os.path.exists('tmp')):
                os.mkdir('tmp')

            pref = self.tab_prefix_ps_check(True)
            fn = 'tmp/{}.{}.pkl'.format(pref, str(rank).zfill(3))

            if os.path.exists(fn) and (not clobber):

                # Should delete if clobber == True?

                if rank == 0:
                    print("Checkpoints for this model found in tmp/.")
                    print("Re-run with clobber=True to overwrite.")

                f = open(fn, 'rb')
                while True:

                    try:
                        tmp = pickle.load(f)
                    except EOFError:
                        break

                    _z.append(tmp[0])
                    _ps.append(tmp[1])
                    _cf.append(tmp[2])

                f.close()

                if _z != []:
                    print("Processor {} loaded checkpoints for z={}".format(rank, _z))

            elif os.path.exists(fn):
                os.remove(fn)

        # Must collect checkpoints so we don't re-run something another
        # processor did!
        if size > 1 and _z != []:
            _zdone = MPI.COMM_WORLD.reduce(_z, root=0)
            zdone = MPI.COMM_WORLD.bcast(_zdone, root=0)
            _zdone_by = MPI.COMM_WORLD.reduce([rank] * len(_z), root=0)
            zdone_by = MPI.COMM_WORLD.bcast(_zdone_by, root=0)
        else:
            zdone = []
            zdone_by = []

        # Figure out what redshift still need to be done by somebody
        assignments = []
        for k, z in enumerate(self.tab_z):
            if z in zdone:
                continue

            assignments.append(z)

        # Split up the work among processors
        my_assignments = []
        for k, z in enumerate(assignments):
            if k % size != rank:
                continue

            my_assignments.append(z)

        if size > 1:
            if len(assignments) % size != 0:
                print("WARNING: Uneven load: {} redshifts and {} processors!".format(len(assignments), size))

        tab_ps_mm = np.zeros((len(self.tab_z), len(self.tab_k)))
        tab_cf_mm = np.zeros((len(self.tab_z), len(self.tab_R)))
        for i, z in enumerate(self.tab_z):

            # Done but not by me!
            if (z in zdone) and (z not in _z):
                continue

            if z not in my_assignments:
                continue

            ##
            # Calculate from scratch
            ##
            print("Processor {} generating z={} PS...".format(rank, z))

            # Must interpolate back to fine grid (uniformly sampled
            # real-space scales) to do FFT and obtain correlation function
            tab_ps_mm[i] = self.get_ps_mm(z, self.tab_k)

            pb.update(i)

            print("Processor {} generating z={} CF...".format(rank, z))

            # Compute correlation function at native resolution to save time
            # later.
            _R_, tab_cf_mm[i] = self.get_cf_mm(z)

            pb.update(i)

            if not checkpoint:
                continue

            with open(fn, 'ab') as f:
                pickle.dump((z, tab_ps_mm[i], tab_cf_mm[i]), f)
                #print("Processor {} wrote checkpoint for z={}".format(rank, z))

        pb.finish()

        # Grab checkpoints before writing to disk
        for i, z in enumerate(self.tab_z):

            # Done but not by me! If not for this, Allreduce would sum
            # solutions from different processors.
            if (z in zdone) and (z not in _z):
                continue

            # Two processors did the same redshift (backward compatibility)
            if zdone.count(z) > 1:
                done_by = []
                for ii, zz in enumerate(zdone):
                    if zz != z:
                        continue
                    done_by.append(zdone_by[ii])

                if rank != done_by[0]:
                    continue

            ##
            # Load checkpoint, if one exists.
            ##
            if z in _z:

                j = _z.index(z)
                tab_ps_mm[i] = _ps[j]
                tab_cf_mm[i] = _cf[j]


        # Collect results!
        if size > 1:
            tmp1 = np.zeros_like(tab_ps_mm)
            nothing = MPI.COMM_WORLD.Allreduce(tab_ps_mm, tmp1)
            self.tab_ps_mm = tmp1

            tmp2 = np.zeros_like(tab_cf_mm)
            nothing = MPI.COMM_WORLD.Allreduce(tab_cf_mm, tmp2)
            self.tab_cf_mm = tmp2

        else:
            self.tab_ps_mm = tab_ps_mm
            self.tab_cf_mm = tab_cf_mm

        # Done!

    def generate_ps(self, fn=None, clobber=True, destination=None, format='hdf5',
        checkpoint=True, **ftkwargs):
        """
        Save matter power spectrum table to HDF5 or binary (via pickle).

        Parameters
        ----------
        fn : str (optional)
            Name of file to save results to. If None, will use
            self.tab_prefix_ps and value of format parameter to make one up.
        clobber : bool
            Overwrite pre-existing files of the same name?
        destination : str
            Path to directory (other than CWD) to save table.
        format : str
            Format of output. Can be 'hdf5' or 'pkl'

        """

        if destination is None:
            destination = '.'

        # Determine filename
        if fn is None:
            fn = '%s/%s.%s' % (destination, self.tab_prefix_ps(True), format)
        else:
            if format not in fn:
                print("Suffix of provided filename does not match chosen format.")
                print("Will go with format indicated by filename suffix.")

        if os.path.exists(fn):
            if clobber:
                os.system('rm -f %s' % fn)
            else:
                raise IOError('File %s exists! Set clobber=True or remove manually.' % fn)

        # Do this first! (Otherwise parallel runs will be garbage)
        self.TabulatePS(clobber=clobber, checkpoint=checkpoint, **ftkwargs)

        if rank > 0:
            return

        self._write_ps(fn, clobber, format)

    def _write_ps(self, fn, clobber, format=format):

        try:
            import hmf
            hmf_v = hmf.__version__
        except AttributeError:
            hmf_v = 'unknown'

        if os.path.exists(fn):
            if clobber:
                os.system('rm -f %s' % fn)
            else:
                raise IOError('File %s exists! Set clobber=True or remove manually.' % fn)

        if format == 'hdf5':
            f = h5py.File(fn, 'w')
            f.create_dataset('tab_z_ps', data=self.tab_z)
            f.create_dataset('tab_R', data=self.tab_R)
            f.create_dataset('tab_k', data=self.tab_k)
            f.create_dataset('tab_ps_mm', data=self.tab_ps_mm)
            f.create_dataset('tab_cf_mm', data=self.tab_cf_mm)

            f.close()
        # Otherwise, pickle it!
        else:
            f = open(fn, 'wb')
            pickle.dump(self.tab_z, f)
            pickle.dump(self.tab_R, f)
            pickle.dump(self.tab_k, f)
            pickle.dump(self.tab_ps_mm, f)
            pickle.dump(self.tab_cf_mm, f)
            pickle.dump(dict(('hmf-version', hmf_v)))
            f.close()

        print('Wrote %s.' % fn)
        return

    def generate_halo_prof(self, format='hdf5', clobber=False, checkpoint=True,
        destination=None, **kwargs):
        """
        Generate a lookup table for Fourier-tranformed halo profiles.
        """

        assert format == 'hdf5'

        if destination is None:
            destination = '.'

        fn = f'{destination}/{self.tab_prefix_prof()}.{format}'

        if rank == 0:
            print(f"# Will save to {fn}.")

        shape = (self.tab_z.size, self.tab_M.size, self.tab_k.size)
        self._tab_u_nfw = np.zeros(shape)

        if self._tab_u_nfw.nbytes / 1e9 > 8:
            print(f"WARNING: Size of profile table projected to be >8 GB!")

        pb = ProgressBar(len(self.tab_z), 'u(z|k,M)', use=rank==0)
        pb.start()

        MM, kk = np.meshgrid(self.tab_M, self.tab_k, indexing='ij')

        for i, z in enumerate(self.tab_z):
            if i % size != rank:
                continue

            self._tab_u_nfw[i,:,:] = self.get_u_nfw(z, MM, kk)
            pb.update(i)

        pb.finish()

        if size > 1:

            tmp = np.zeros(shape)
            nothing = MPI.COMM_WORLD.Allreduce(self._tab_u_nfw, tmp)
            self._tab_u_nfw = tmp

            # So only root processor writes to disk
            if rank > 0:
                return

        with h5py.File(fn, 'w') as f:
            f.create_dataset('tab_u_nfw', data=self._tab_u_nfw)
            f.create_dataset('tab_k', data=self.tab_k)
            f.create_dataset('tab_M', data=self.tab_M)
            f.create_dataset('tab_z', data=self.tab_z)

        if rank == 0:
            print(f"# Wrote {fn}.")

        return

    def generate_halo_surface_dens(self, format='hdf5', clobber=False,
        checkpoint=True, destination=None, **kwargs):
        """
        Generate a lookup table for halo surface density.
        """

        assert format == 'hdf5'

        if destination is None:
            destination = '.'

        fn = f'{destination}/{self.tab_prefix_surf()}.{format}'
        if rank == 0:
            print(f"# Will save to {fn}.")

        # Hard-coded for now, change this.
        Rmi, Rma = -3, 1
        dlogR = 0.25
        R = 10**np.arange(Rmi, Rma+dlogR, dlogR)

        shape = (self.tab_z.size, self.tab_M.size, R.size)
        self._tab_sigma_nfw = np.zeros(shape)
        if self._tab_sigma_nfw.nbytes / 1e9 > 8:
            print(f"WARNING: Size of profile table projected to be >8 GB!")

        pb = ProgressBar(len(self.tab_z), 'Sigma(z|M,R)', use=rank==0)
        pb.start()

        for i, z in enumerate(self.tab_z):
            if i % size != rank:
                continue

            pb.update(i)

            model_nfw = lambda MM, rr: self.get_rho_nfw(z, Mh=MM, r=rr,
                truncate=False)

            # Tabulate surface density as a function of displacement
            # and halo mass

            dlogm = self.pf['halo_dlogM']
            Sall = np.zeros((self.tab_M.size, R.size))

            for ii, _M_ in enumerate(self.tab_M):
                rho = lambda rr: model_nfw(_M_, rr)
                Sigma = lambda R: 2 * \
                    quad(lambda rr: rr * rho(rr) / np.sqrt(rr**2 - R**2),
                    R, np.inf)[0]
                for jj, _R_ in enumerate(R):
                    self._tab_sigma_nfw[i,ii,jj] = Sigma(_R_)

        pb.finish()

        if size > 1:

            tmp = np.zeros(shape)
            nothing = MPI.COMM_WORLD.Allreduce(self._tab_sigma_nfw, tmp)
            self._tab_sigma_nfw = tmp

            # So only root processor writes to disk
            if rank > 0:
                return

        with h5py.File(fn, 'w') as f:
            f.create_dataset('tab_Sigma_nfw', data=self._tab_sigma_nfw)
            f.create_dataset('tab_R', data=R)
            f.create_dataset('tab_M', data=self.tab_M)
            f.create_dataset('tab_z', data=self.tab_z)

        if rank == 0:
            print(f"# Wrote {fn}.")
