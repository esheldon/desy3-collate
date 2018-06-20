"""
collate the metacal files, renaming columns to follow
a standard.

also calculate some new columns
"""
# note I set the defaults for integers to -9999 rather than -1

from __future__ import print_function
import os
import numpy as np
import fitsio
import glob
import esutil as eu
import time

from . import util
from .util import get_namer

# some data and data structures
from .util import (
    DEFAULT,
    PDEFAULT,
    MCAL_SN,
    MCAL_NAME_MAP,
    MCAL_SHEAR_NAME_MAP,
    MCAL_R_DT,
    MCAL_SHEAR_TYPES,
    dgamma,
)

class Collator(object):
    """
    collate the tiles, renaming fields to a common standard.
    """
    def __init__(self, metacal_dir, output_file):
        self._metacal_dir=metacal_dir
        self._output_file=output_file
        self._load_flist()
    
    def go(self):
        """
        go by tile and write the data
        """
        nfiles=len(self._flist)

        print("writing to:",self._output_file)

        first=True
        with fitsio.FITS(self._output_file,'rw',clobber=True) as fits:
            for i,fname in enumerate(self._flist):

                print("%d/%d %s" % (i+1,nfiles,fname))

                print("    reading")
                data = fitsio.read(fname)

                print("    extracting output")

                odata = self._extract_output(data)

                if first:
                    first=False
                    fits.write(odata,extname=util.EXTNAME)
                else:
                    fits[-1].append(odata)

                #if i==4:
                #    return
        print("output is in:",self._output_file)

    def _extract_output(self, data): 
        """
        extract the outputs, following the name mapping.

        Calculate any extra fields
        """

        w,=np.where(data['flags']==0)

        output = self._get_collated_struct(data)

        nbands = util.get_nbands_from_pars(data['mcal_pars'][0])

        nmap = util.get_mcal_name_map(nbands)
        for name,nm in nmap.iteritems():
            copy_data = self._get_copy_data(data, nm['name'], nm['el'])
            output[name] = copy_data

        nmap = util.get_mcal_shear_name_map(nbands)

        for t_oname,nm in nmap.iteritems():
            for type in MCAL_SHEAR_TYPES:
                if self._should_skip_mcal_type(nm, type):
                    continue

                namer = get_namer(type=type)

                oname = namer( t_oname )
                dname = namer( nm['name'] )

                copy_data = self._get_copy_data(data, dname, nm['el'])

                # we need to sqrt the covariance elements in some
                # cases
                if self._is_calculated_err(oname, dname):
                    self._convert_err_from_cov(copy_data)

                output[oname] = copy_data
            

        if w.size > 0:
            output['R11'][w] = (output['e1_1p'][w] - output['e1_1m'][w])/dgamma
            output['R12'][w] = (output['e1_2p'][w] - output['e1_2m'][w])/dgamma
            output['R21'][w] = (output['e2_1p'][w] - output['e2_1m'][w])/dgamma
            output['R22'][w] = (output['e2_2p'][w] - output['e2_2m'][w])/dgamma

        self._fix_data(output)
        return output

    def _get_copy_data(self, data, name, el):
        """
        deal with non-scalar data
        """
        if el is None:
            copy_data = data[name]
        elif isinstance(el,tuple):
            el1,el2 = el
            copy_data = data[name][:,el1,el2]
        else:
            copy_data = data[name][:,el]

        return copy_data

    def _fix_data(self, data):
        """
        fix columns that are known to have NaN etc.
        """
        for name in ['psf_e1','psf_e2','psf_size']:
            w,=np.where(np.isfinite(data[name])==False)
            if w.size > 0:
                #print("        fixing %d NaN in %s" % (w.size,name))
                data[name][w] = -9999

    def _load_flist(self):
        """
        get the list of files
        """
        self._flist = util.load_mcal_flist(self._metacal_dir)

    def _get_collated_struct(self, data):
        nbands = util.get_nbands_from_pars(data['mcal_pars'][0])

        dt = self._get_collated_dtype(nbands)
        st = np.zeros(data.size, dtype=dt)

        util.set_defaults(st)
        return st

    def _get_collated_dtype(self, nbands):
        """
        build up the output data type from the
        name map dict
        """
        dt = []

        nmap = util.get_mcal_name_map(nbands)
        for name,nm in nmap.iteritems():
            dt.append( (name, nm['dt']) )

        nmap = util.get_mcal_shear_name_map(nbands)
        for name,nm in nmap.iteritems():
            for type in MCAL_SHEAR_TYPES:
                if self._should_skip_mcal_type(nm, type):
                    continue
                namer = get_namer(type=type)

                mcal_name = namer(name)
                dt.append( (mcal_name, nm['dt']) )

        dt += MCAL_R_DT
        return dt

    def _is_calculated_err(self, oname, dname):
        if '_err' in oname and '_cov' in dname:
            return True
        else:
            return False

    def _convert_err_from_cov(self, cov):
        w,=np.where(cov > 0)
        if w.size > 0:
            cov[w] = np.sqrt(cov[w])

    def _should_skip_mcal_type(self, nmap, type):
        skip=False
        if 'mcal_types' in nmap:
            if type not in nmap['mcal_types']:
                #print("    skipping type:",type)
                skip=True 

        return skip
