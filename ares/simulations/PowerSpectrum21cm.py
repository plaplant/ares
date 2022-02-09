import os
import copy
import pickle
import numpy as np
from types import FunctionType
from ..static import Fluctuations
from .Global21cm import Global21cm
from ..physics.HaloModel import HaloModel
from ..util import ParameterFile, ProgressBar
#from ..analysis.BlobFactory import BlobFactory
from ..physics.Constants import cm_per_mpc, c, s_per_yr
from ..analysis.PowerSpectrum import PowerSpectrum as AnalyzePS

#
#try:
#    import dill as pickle
#except ImportError:
#    import pickle

defaults = \
{
 'load_ics': True,
}

class PowerSpectrum21cm(AnalyzePS): # pragma: no cover
    def __init__(self, **kwargs):
        """ Set up a power spectrum calculation. """

        # See if this is a tanh model calculation
        #is_phenom = self._check_if_phenom(**kwargs)

        #kwargs.update(defaults)
        #if 'problem_type' not in kwargs:
        #    kwargs['problem_type'] = 101

        self.kwargs = kwargs

    @property
    def mean_history(self):
        if not hasattr(self, '_mean_history'):
            self.gs.run()
            self._mean_history = self.gs.history

        return self._mean_history

    @mean_history.setter
    def mean_history(self, value):
        self._mean_history = value

    @property
    def pops(self):
        return self.gs.medium.field.pops

    @property
    def grid(self):
        return self.gs.medium.field.grid

    @property
    def hydr(self):
        return self.grid.hydr

    @property
    def pf(self):
        if not hasattr(self, '_pf'):
            self._pf = ParameterFile(**self.kwargs)
        return self._pf

    @pf.setter
    def pf(self, value):
        self._pf = value

    @property
    def gs(self):
        if not hasattr(self, '_gs'):
            self._gs = Global21cm(**self.kwargs)
        return self._gs

    @gs.setter
    def gs(self, value):
        """ Set global 21cm instance by hand. """
        self._gs = value

    @property
    def field(self):
        if not hasattr(self, '_field'):
            self._field = Fluctuations(**self.kwargs)
        return self._field

    @property
    def halos(self):
        if not hasattr(self, '_halos'):
            self._halos = self.pops[0].halos
        return self._halos

    @property
    def tab_z(self):
        if not hasattr(self, '_tab_z'):
            self._tab_z = np.array(np.sort(self.pf['ps_output_z'])[-1::-1],
                dtype=np.float64)
        return self._tab_z

    @tab_z.setter
    def tab_z(self, value):
        if type(value) == np.ndarray:
            self._tab_z = value
        else:
            self._tab_z = np.array(value)

    def run(self, z=None, k=None):
        """
        Run a simulation, compute power spectrum at each redshift.

        Returns
        -------
        Nothing: sets `history` attribute.

        """

        if z is not None:
            self.tab_z = z
        if k is not None:
            self.tab_k = k

        N = self.tab_z.size

        pb = self.pb = ProgressBar(N, use=self.pf['progress_bar'],
            name='ps-21cm')

        all_ps = []
        for i, (z, data) in enumerate(self.step()):

            # Do stuff
            all_ps.append(data.copy())

            if i == 0:
                keys = data.keys()

            if not pb.has_pb:
                pb.start()

            pb.update(i)

        pb.finish()

        self.all_ps = all_ps

        hist = {}
        for key in keys:

            is2d_k = key.startswith('ps')
            is2d_R = key.startswith('jp') or key.startswith('ev') \
                  or key.startswith('cf')
            is2d_B = (key in ['n_i', 'm_i', 'r_i', 'delta_B'])

            if is2d_k:
                tmp = np.zeros((len(self.tab_z), len(self.tab_k)))
            elif is2d_R:
                tmp = np.zeros((len(self.tab_z), len(self.tab_R)))
            elif is2d_B:
                tmp = np.zeros((len(self.tab_z), len(all_ps[0]['r_i'])))
            else:
                tmp = np.zeros_like(self.tab_z)

            for i, z in enumerate(self.tab_z):
                if key not in all_ps[i].keys():
                    continue

                tmp[i] = all_ps[i][key]

            hist[key] = tmp.copy()

        self.history = hist
        self.history['z'] = self.tab_z
        self.history['k'] = self.tab_k
        self.history['R'] = self.tab_R

    @property
    def tab_k(self):
        """
        Wavenumbers to output power spectra.

        .. note :: Can be far more crude than native resolution of
            matter power spectrum.

        """

        if not hasattr(self, '_tab_k'):
            if self.pf['ps_output_k'] is not None:
                self._tab_k = self.pf['ps_output_k']
            else:
                lnk1 = self.pf['ps_output_lnkmin']
                lnk2 = self.pf['ps_output_lnkmax']
                dlnk = self.pf['ps_output_dlnk']
                self._tab_k = np.exp(np.arange(lnk1, lnk2+dlnk, dlnk))

        return self._tab_k

    @tab_k.setter
    def tab_k(self, value):
        if type(value) == np.ndarray:
            self._tab_k = value
        else:
            self._tab_k = np.array(value)

    @property
    def tab_R(self):
        """
        Scales on which to compute correlation functions.

        .. note :: Can be more crude than native resolution of matter
            power spectrum, however, unlike `self.tab_k`, the resolution of
            this quantity matters when converting back to power spectra,
            since that operation requires an integral over R.

        """
        if not hasattr(self, '_R'):
            if self.pf['ps_output_R'] is not None:
                self._R = self.pf['ps_output_R']
            else:
                lnR1 = self.pf['ps_output_lnRmin']
                lnR2 = self.pf['ps_output_lnRmax']
                dlnR = self.pf['ps_output_dlnR']
                #lnR = np.log(self.halos.tab_R)

                self._R = np.exp(np.arange(lnR1, lnR2+dlnR, dlnR))

        return self._R

    @property
    def tab_Mmin(self):
        if not hasattr(self, '_tab_Mmin'):
            self._tab_Mmin = np.ones_like(self.halos.tab_z) * np.inf
            for j, pop in enumerate(self.pops):
                self._tab_Mmin = np.minimum(self._tab_Mmin, pop._tab_Mmin)

        return self._tab_Mmin

    @property
    def tab_zeta(self):
        pass

    def step(self):
        """
        Generator for the power spectrum.
        """

        transform_kwargs = dict(split_by_scale=self.pf['ps_split_transform'],
                epsrel=self.pf['ps_fht_rtol'],
                epsabs=self.pf['ps_fht_atol'])

        # Set a few things before we get moving.
        self.field.tab_Mmin = self.tab_Mmin

        for i, z in enumerate(self.tab_z):

            data = {}

            ##
            # First, loop over populations and determine total
            # UV and X-ray outputs.
            ##

            # Prepare for the general case of Mh-dependent things
            Nion = np.zeros_like(self.halos.tab_M)
            Nlya = np.zeros_like(self.halos.tab_M)
            fXcX = np.zeros_like(self.halos.tab_M)
            zeta_ion = zeta = np.zeros_like(self.halos.tab_M)
            zeta_lya = np.zeros_like(self.halos.tab_M)
            zeta_X = np.zeros_like(self.halos.tab_M)
            #Tpro = None
            for j, pop in enumerate(self.pops):
                pop_zeta = pop.get_zeta_ion(z=z)
                zeta += pop_zeta

            ##
            # Make scalar if it's a simple model
            ##
            if np.all(np.diff(zeta) == 0):
                zeta = zeta[0]
            if np.all(np.diff(zeta_X) == 0):
                zeta_X = zeta_X[0]
            if np.all(np.diff(zeta_lya) == 0):
                zeta_lya = zeta_lya[0]

            self.field.zeta = zeta
            self.field.zeta_X = zeta_X

            self.zeta = zeta

            ##
            # Figure out scaling from ionized regions to heated regions.
            # Right now, only constant (relative) scaling is allowed.
            ##
            asize = self.pf['bubble_shell_asize_zone_0']
            if self.pf['ps_include_temp'] and asize is not None:

                self.field.is_Rs_const = False

                if type(asize) is FunctionType:
                    R_s = lambda R, z: R + asize(z)
                else:
                    R_s = lambda R, z: R + asize

            elif self.pf['ps_include_temp'] and self.pf['ps_include_ion']:
                fvol = self.pf["bubble_shell_rvol_zone_0"]
                frad = self.pf['bubble_shell_rsize_zone_0']

                assert (fvol is not None) + (frad is not None) <= 1

                if fvol is not None:
                    assert frad is None

                    # Assume independent variable is redshift for now.
                    if type(fvol) is FunctionType:
                        frad = lambda z: (1. + fvol(z))**(1./3.) - 1.
                        self.field.is_Rs_const = False
                    else:
                        frad = lambda z: (1. + fvol)**(1./3.) - 1.

                elif frad is not None:
                    if type(frad) is FunctionType:
                        self.field.is_Rs_const = False
                    else:
                        frad = lambda z: frad
                else:
                    # If R_s = R_s(z), must re-compute overlap volumes on each
                    # step. Should set attribute if this is the case.
                    raise NotImplemented('help')

                R_s = lambda R, z: R * (1. + frad(z))


            else:
                R_s = lambda R, z: None
                Th = None

            # Must be constant, for now.
            Th = self.pf["bubble_shell_ktemp_zone_0"]

            self.R_s = R_s
            self.Th = Th


            ##
            # First: some global quantities we'll need
            ##
            Tcmb = self.cosm.TCMB(z)
            Tk = np.interp(z, self.mean_history['z'][-1::-1],
                self.mean_history['igm_Tk'][-1::-1])
            Ts = np.interp(z, self.mean_history['z'][-1::-1],
                self.mean_history['Ts'][-1::-1])
            Ja = np.interp(z, self.mean_history['z'][-1::-1],
                self.mean_history['Ja'][-1::-1])
            xHII, ne = [0] * 2

            xa = self.hydr.RadiativeCouplingCoefficient(z, Ja, Tk)
            xc = self.hydr.CollisionalCouplingCoefficient(z, Tk)
            xt = xa + xc

            # Won't be terribly meaningful if temp fluctuations are off.
            C = self.field.TempToContrast(z, Th=Th, Tk=Tk, Ts=Ts, Ja=Ja)
            data['c'] = C
            data['Ts'] = Ts
            data['Tk'] = Tk
            data['xa'] = xa
            data['Ja'] = Ja

            # Assumes strong coupling. Mapping between temperature
            # fluctuations and contrast fluctuations.
            #Ts = Tk

            # Add beta factors to dictionary
            for f1 in ['x', 'd', 'a']:
                func = self.hydr.__getattribute__('beta_%s' % f1)
                data['beta_%s' % f1] = func(z, Tk, xHII, ne, Ja)

            Qi_gs = np.interp(z, self.gs.history['z'][-1::-1],
                self.gs.history['cgm_h_2'][-1::-1])

            # Ionization fluctuations
            if self.pf['ps_include_ion']:

                Ri, Mi, Ni = self.field.BubbleSizeDistribution(z, ion=True)

                data['n_i'] = Ni
                data['m_i'] = Mi
                data['r_i'] = Ri
                data['delta_B'] = self.field._B(z, ion=True)
            else:
                Ri = Mi = Ni = None

            Qi = self.field.MeanIonizedFraction(z)

            Qi_bff = self.field.BubbleFillingFactor(z)

            xibar = Qi_gs

            #print(z, Qi_bff, Qi, xibar, Qi_bff / Qi)

            if self.pf['ps_include_temp']:
                # R_s=R_s(Ri,z)
                Qh = self.field.MeanIonizedFraction(z, ion=False)
                data['Qh'] = Qh
            else:
                data['Qh'] = Qh = 0.0

            # Interpolate global signal onto new (coarser) redshift grid.
            dTb_ps = np.interp(z, self.gs.history['z'][-1::-1],
                self.gs.history['dTb'][-1::-1])

            xavg_gs = np.interp(z, self.gs.history['z'][-1::-1],
                self.gs.history['xavg'][-1::-1])

            data['dTb'] = dTb_ps

            #data['dTb_bulk'] = np.interp(z, self.gs.history['z'][-1::-1],
            #    self.gs.history['dTb_bulk'][-1::-1])


            ##
            # Correct for fraction of ionized and heated volumes
            # and densities!
            ##
            if self.pf['ps_include_temp']:
                data['dTb_vcorr'] = None#(1 - Qh - Qi) * data['dTb_bulk'] \
                    #+ Qh * self.hydr.dTb(z, 0.0, Th)
            else:
                data['dTb_vcorr'] = None#data['dTb_bulk'] * (1. - Qi)

            if self.pf['ps_include_xcorr_ion_rho']:
                pass
            if self.pf['ps_include_xcorr_ion_hot']:
                pass

            # Just for now
            data['dTb0'] = data['dTb']
            data['dTb0_2'] = data['dTb0_1'] = data['dTb_vcorr']

            #if self.pf['include_ion_fl']:
            #    if self.pf['ps_rescale_Qion']:
            #        xibar = min(np.interp(z, self.pops[0].halos.z,
            #            self.pops[0].halos.fcoll_Tmin) * zeta, 1.)
            #        Qi = xibar
            #
            #        xibar = np.interp(z, self.mean_history['z'][-1::-1],
            #            self.mean_history['cgm_h_2'][-1::-1])
            #
            #    else:
            #        Qi = self.field.BubbleFillingFactor(z, zeta)
            #        xibar = 1. - np.exp(-Qi)
            #else:
            #    Qi = 0.



            #if self.pf['ps_force_QHII_gs'] or self.pf['ps_force_QHII_fcoll']:
            #    rescale_Q = True
            #else:
            #    rescale_Q = False

            #Qi = np.mean([QHII_gs, self.field.BubbleFillingFactor(z, zeta)])

            #xibar = np.interp(z, self.mean_history['z'][-1::-1],
            #    self.mean_history['cgm_h_2'][-1::-1])

            # Avoid divide by zeros when reionization is over
            if Qi == 1:
                Tbar = 0.0
            else:
                Tbar = data['dTb0_2']

            xbar = 1. - xibar
            data['Qi'] = Qi
            data['xibar'] = xibar
            data['dTb0'] = Tbar
            #data['dTb_bulk'] = dTb_ps / (1. - xavg_gs)

            ##
            # 21-cm fluctuations
            ##
            data['cf_21'] = self.field.get_cf(z,
                R=self.tab_R, term='21', R_s=R_s(Ri,z), Ts=Ts, Th=Th,
                Tk=Tk, Ja=Ja, k=self.tab_k)

            # Always compute the 21-cm power spectrum. Individual power
            # spectra can be saved by setting ps_save_components=True.
            data['ps_21'] = self.field.get_ps_from_cf(self.tab_k,
                data['cf_21'], R=self.tab_R, **transform_kwargs)

            # Always save the matter correlation function.
            data['cf_dd'] = self.field.get_cf(z, term='dd', R=self.tab_R)

            # Grab the ionization CF
            data['cf_bb'] = self.field.get_cf(z, term='ii', R=self.tab_R)

            if self.pf['ps_output_components']:
                data['ps_dd'] = self.field.get_ps_from_cf(self.tab_k,
                    data['cf_dd'], R=self.tab_R, **transform_kwargs)
                data['ps_bb'] = self.field.get_ps_from_cf(self.tab_k,
                    data['cf_bb'], R=self.tab_R, **transform_kwargs)





            yield z, data

    def save(self, prefix, suffix='pkl', clobber=False, fields=None):
        """
        Save results of calculation. Pickle parameter file dict.

        Notes
        -----
        1) will save files as prefix.history.suffix and prefix.parameters.pkl.
        2) ASCII files will fail if simulation had multiple populations.

        Parameters
        ----------
        prefix : str
            Prefix of save filename
        suffix : str
            Suffix of save filename. Can be hdf5 (or h5) or pkl.
            Anything else will be assumed to be ASCII format (e.g., .txt).
        clobber : bool
            Overwrite pre-existing files of same name?

        """

        self.gs.save(prefix, clobber=clobber, fields=fields)

        fn = '%s.fluctuations.%s' % (prefix, suffix)

        if os.path.exists(fn):
            if clobber:
                os.remove(fn)
            else:
                raise IOError('%s exists! Set clobber=True to overwrite.' % fn)

        if suffix == 'pkl':
            f = open(fn, 'wb')
            pickle.dump(self.history._data, f)
            f.close()

            try:
                f = open('%s.blobs.%s' % (prefix, suffix), 'wb')
                pickle.dump(self.blobs, f)
                f.close()

                if self.pf['verbose']:
                    print('Wrote {}.blobs.{}'.format(prefix, suffix))
            except AttributeError:
                print('Error writing {}.blobs.{}'.format(prefix, suffix))

        elif suffix in ['hdf5', 'h5']:
            import h5py

            f = h5py.File(fn, 'w')
            for key in self.history:
                if fields is not None:
                    if key not in fields:
                        continue
                f.create_dataset(key, data=np.array(self.history[key]))
            f.close()

        # ASCII format
        else:
            f = open(fn, 'w')
            print >> f, "#",

            for key in self.history:
                if fields is not None:
                    if key not in fields:
                        continue
                print >> f, '%-18s' % key,

            print >> f, ''

            # Now, the data
            for i in range(len(self.history[key])):
                s = ''

                for key in self.history:
                    if fields is not None:
                        if key not in fields:
                            continue

                    s += '%-20.8e' % (self.history[key][i])

                if not s.strip():
                    continue

                print >> f, s

            f.close()

        if self.pf['verbose']:
            print('Wrote {}.fluctuations.{}'.format(prefix, suffix))

        #write_pf = True
        #if os.path.exists('%s.parameters.pkl' % prefix):
        #    if clobber:
        #        os.remove('%s.parameters.pkl' % prefix)
        #    else:
        #        write_pf = False
        #        print 'WARNING: %s.parameters.pkl exists! Set clobber=True to overwrite.' % prefix

        #if write_pf:
        #
        #    #pf = {}
        #    #for key in self.pf:
        #    #    if key in self.carryover_kwargs():
        #    #        continue
        #    #    pf[key] = self.pf[key]
        #
        #    if 'revision' not in self.pf:
        #        self.pf['revision'] = get_hg_rev()
        #
        #    # Save parameter file
        #    f = open('%s.parameters.pkl' % prefix, 'wb')
        #    pickle.dump(self.pf, f, -1)
        #    f.close()
        #
        #    if self.pf['verbose']:
        #        print 'Wrote %s.parameters.pkl' % prefix
        #
