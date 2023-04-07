"""

Misc.py

Author: Jordan Mirocha
Affiliation: University of Colorado at Boulder
Created on: Sun Oct 19 19:50:31 MDT 2014

Description:

"""
import os
import subprocess
import numpy as np
from ..data import ARES
from .Stats import bin_e2c

numeric_types = [int, float, np.int64, np.float64]

def get_cmd_line_kwargs(argv):

    cmd_line_kwargs = {}

    for arg in argv[1:]:
        try:
            pre, post = arg.split('=')
        except ValueError:
            # To deal with parameter values that have an '=' in them.
            pre = arg[0:arg.find('=')]
            post = arg[arg.find('=')+1:]

        # Need to do some type-casting
        if post.isdigit():
            cmd_line_kwargs[pre] = int(post)
        elif post.isalpha():
            if post == 'None':
                cmd_line_kwargs[pre] = None
            elif post in ['True', 'False']:
                cmd_line_kwargs[pre] = True if post == 'True' else False
            else:
                cmd_line_kwargs[pre] = str(post)
        elif post[0] == '[':
            vals = post[1:-1].split(',')
            cmd_line_kwargs[pre] = np.array([float(val) for val in vals])
        else:
            try:
                cmd_line_kwargs[pre] = float(post)
            except ValueError:
                # strings with underscores will return False from isalpha
                cmd_line_kwargs[pre] = str(post)

    return cmd_line_kwargs

def get_hash(repo_path=ARES, repo_env=None):
    """
    Return the unique git hash associated with the HEAD of some repository.

    This intended to be used to save the current version of some code you're
    using with any output files to help with debugging later. For example,
    I have some output files for a calculation I've done that change over time,
    indicating a problem/development, and I want to know precisely when that
    change happened. In practice, I usually save the output of this function
    as metadata in an hdf5 file (e.g., as an attribute of some dataset or
    as a dataset of its own).

    Parameters
    ----------
    repo_path : str, None
        Absolute path to root directory of repo of interest (where .git lives)
    repo_env : str, None
        Name of environment variable that points to repo (again, the root
        directory where .git lives).

    Returns
    -------
    A string containing the unique hash of the current HEAD of git repo.

    Known Flaws
    -----------
    If you run your code with uncommitted changes, then this hash may not help
    you find a bug, as the bug could have been in your uncommitted changes.
    There's not really a good solution to this, other than to always run your
    code with a 'clean' install!

    """

    assert (repo_path is not None) or (repo_env is not None), \
        "Must supply path to git repo or environment variable that points to it."

    try:
        cwd = os.getcwd()

        if repo_env is not None:
            PATH = os.environ.get(repo_env)
        else:
            PATH = repo_path

        os.chdir(PATH)

        # git rev-parse HEAD
        pipe = subprocess.Popen(["git", "rev-parse", "HEAD"],
            stdout=subprocess.PIPE)

        # Move back to where we were
        os.chdir(cwd)
    except Exception as err:
        print("Failure to obtain hash due to following error: {}".format(err))
        return 'unknown'

    return pipe.stdout.read().strip()

def num_freq_bins(Nx, zi=40, zf=10, Emin=2e2, Emax=3e4):
    """
    Compute number of frequency bins required for given log-x grid.

    Defining the variable x = 1 + z, and setting up a grid in log-x containing
    Nx elements, compute the number of frequency bins required for a 1:1
    mapping between redshift and frequency.

    """
    x = np.logspace(np.log10(1.+zf), np.log10(1.+zi), Nx)
    R = x[1] / x[0]

    # Create mapping to frequency space
    Etmp = 1. * Emin
    n = 1
    while Etmp < Emax:
        Etmp = Emin * R**(n - 1)
        n += 1

    # Subtract 2: 1 because we overshoot Emax in while loop, another because
    # n is index-1-based (?)

    return n-2

def get_attribute(s, ob):
    """
    Break apart a string `s` and recursively fetch attributes from object `ob`.
    """
    spart = s.partition('.')

    f = ob
    for part in spart:
        if part == '.':
            continue

        f = f.__getattribute__(part)

    return f

def split_by_sign(x, y):
    """
    Split apart an array into its positive and negative chunks.
    """

    splitter = np.diff(np.sign(y))

    if np.all(splitter == 0):
        ych = [y]
        xch = [x]
    else:
        splits = np.atleast_1d(np.argwhere(splitter != 0).squeeze()) + 1
        ych = np.split(y, splits)
        xch = np.split(x, splits)

    return xch, ych

def get_field_from_catalog(field, pos, Lbox, dims=512, mesh=None,
    weight_by_field=True, by_volume=True):
    """
    Convert a catalog, i.e., a list of (lum, x, y, z), to luminosity
    (or whatever) on a mesh.

    .. note :: If you're applying some threshold like Mmin, do so
        BEFORE running this routine. In this case, ``catalog`` should
        be a numpy masked array.

    .. note :: If weight_by_field == False, the units of the output
        will just number of halos per voxel, i.e., independent
        of the field.

    Parameters
    ----------
    catalog : np.ndarray
        Should have shape (Ngalaxies, 4)
    dims : int
        Linear dimensions of box to create. Can alternatively provide
        desired grid resolution (in Mpc / h) via ``mesh`` keyword
        argument (see below).
    mesh : int, float
        If supplied, should be the linear dimension of voxels used
        in histogram [Mpc / h].
    weight_by_field : bool
        If True, will weight by field (density or luminosity usually)
    by_volume : bool
        If True, will divide by voxel volume so field has units of
        x / cMpc^3, where x = whatever the field is (e.g., mass, luminosity).
        Otherwise, units will be the same as the input array. We generally
        set this to True for things like halo mass density, and False
        for things like the total ionizing photon output, which we want
        as an absolute photon production rate, not production rate density.

    """

    if mesh is None:
        mesh = Lbox / float(dims)

    xe = np.arange(0, Lbox+mesh, mesh)
    ye = np.arange(0, Lbox+mesh, mesh)
    ze = np.arange(0, Lbox+mesh, mesh)

    _x, _y, _z = pos.T

    data = np.array([_x, _y, _z]).T
    hist, edges = np.histogramdd(data, bins=[xe, ye, ze],
        weights=field if weight_by_field else None, normed=False)

    if by_volume:
        hist /= mesh**3

    return bin_e2c(xe), hist
