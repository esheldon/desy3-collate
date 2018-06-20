"""
test a random set of tiles

The collation code uses an algorithmic name map, so here I intentionally used a
different way to transform the names in order to test the name map.  The
new method could have bugs in it as well

"""

from __future__ import print_function
import os
import numpy
import fitsio
import glob
import esutil as eu
import random

from . import util
from .util import dgamma

class Namer(object):
    """
    wrote from scratch to test Namer in util
    """
    def __init__(self, type=None):
        if type=='noshear':
            type=None

        self._type=type

    def __call__(self, name):
        if self._type is not None:
            n = '%s_%s' % (name,self._type)
        else:
            n=name

        return n

class Tester(object):
    """
    test we copied and calculated things correctly
    """
    def __init__(self, metacal_dir, shear_file, ntest=3, seed=None):
        self._metacal_dir=metacal_dir
        self._shear_file=shear_file
        self._ntest=ntest

        if seed is not None:
            random.seed(seed)

        self._load_flist()
        self._load_ids()
    
    def go(self):
        """
        test each of the random tiles
        """
        nfiles=len(self._flist)

        for i,fname in enumerate(self._flist):

            print("%d/%d %s" % (i+1,nfiles,fname))

            print("    reading")
            data = fitsio.read(fname)

            print("    matching")
            mdata, matches = self._get_matches(data)

            if mdata.size == 0:
                continue

            print("    testing")
            self._test(data[mdata], matches)

    def _test(self, data, matches):
        """
        Intentionally hard coding some names and indices as a test of the
        algorithmic method in the collator
        """
        from numpy import all

        nbands = util.get_nbands_from_pars(data['mcal_pars'][0])
        bands=util.get_bands(nbands)

        for n in ['flags','mask_frac']:
            assert all(data[n] == matches[n])

        for type in util.MCAL_SHEAR_TYPES:
            n=Namer(type=type)
                
            assert all(data[n('mcal_g')][:,0] == matches[n('e1')])
            assert all(data[n('mcal_g')][:,1] == matches[n('e2')])

            assert all(data[n('mcal_T_r')] == matches[n('T')])
            assert all(data[n('mcal_s2n_r')] == matches[n('snr')])

            pn=n('mcal_pars')
            for i,band in enumerate(bands):
                bn = n('flux_%s' % band)
                assert all(data[pn][:,5+i]==matches[bn])

            cn=n('mcal_g_cov')
            assert all(data[cn][:,0,0] == matches[n('covmat_0_0')])
            assert all(data[cn][:,0,1] == matches[n('covmat_0_1')])
            assert all(data[cn][:,1,1] == matches[n('covmat_1_1')])

        w,=numpy.where(data['flags'] == 0)
        if w.size > 0:
            R11 = (data['mcal_g_1p'][w,0] - data['mcal_g_1m'][w,0])/dgamma
            R12 = (data['mcal_g_2p'][w,0] - data['mcal_g_2m'][w,0])/dgamma
            R21 = (data['mcal_g_1p'][w,1] - data['mcal_g_1m'][w,1])/dgamma
            R22 = (data['mcal_g_2p'][w,1] - data['mcal_g_2m'][w,1])/dgamma

            assert all(R11 == matches['R11'][w])
            assert all(R12 == matches['R12'][w])
            assert all(R21 == matches['R21'][w])
            assert all(R22 == matches['R22'][w])

    def _get_matches(self, data):
        mdata, mids = eu.numpy_util.match(
            data['id'],
            self._ids,
        )

        if mdata.size == 0:
            print("    none matched!")
            return mdata, None
        
        matches = fitsio.read(
            self._shear_file,
            rows=mids,
        )

        return mdata, matches


    def _load_ids(self):
        print("loading ids from collated catalog:",
              self._shear_file)
        self._ids = fitsio.read(
            self._shear_file,
            columns='coadd_objects_id',
        )

    def _load_flist(self):
        """
        get the list of files
        """
        print("loading file list")
        pattern = 'DES*blind.fits'
        pattern = os.path.join(self._metacal_dir, pattern)

        self._flist = glob.glob(pattern)
        if len(self._flist) == 0:
            raise RuntimeError("found no files '%s'" % pattern)

        random.shuffle(self._flist)
        self._flist = self._flist[0:self._ntest]

