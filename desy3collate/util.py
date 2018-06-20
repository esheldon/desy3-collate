"""
notes

    reorder to match all to main gold cat
        - combine collated region files first

collation of all regions into one big file
    region column

negation of flags into mcal_select field
"""
from __future__ import print_function
from collections import OrderedDict
import os
import glob
import copy

CHUNKSIZE = 1000000

DEFAULT = -9999
PDEFAULT = 9999
SDEFAULT = 'None'

FLAG_NOT_MEASURED = 2**30

MOF_BANDS=['g','r','i','z']

MOF_BADFLAGS = 2**8 + 2**9
MCAL_SN=0.2
dgamma = 2*0.01

# BPZ cannot deal with negative fluxes
# this corresponds to a mag of 45
MINFLUX = 1.0e-6

EXTNAME='model_fits'

MCAL_FILE_PATTERN='DES*blind.fits'

MCAL_NAME_MAP=OrderedDict([
    ('coadd_objects_id', {'name':'id',         'dt':'i8', 'el':None}),
    ('flags',            {'name':'flags',      'dt':'i4', 'el':None}),
    ('mask_frac',        {'name':'mask_frac',  'dt':'f8', 'el':None}),

    ('psf_e1',           {'name':'psfrec_g',   'dt':'f8', 'el':0}),
    ('psf_e2',           {'name':'psfrec_g',   'dt':'f8', 'el':1}),
    ('psf_size',         {'name':'psfrec_T',   'dt':'f8', 'el':None}),

    ('mcal_psf_e1',      {'name':'mcal_gpsf',  'dt':'f8', 'el':0}),
    ('mcal_psf_e2',      {'name':'mcal_gpsf',  'dt':'f8', 'el':1}),
    ('mcal_psf_size',    {'name':'mcal_Tpsf',  'dt':'f8', 'el':None}),
])

DEC_MIN=-60.0
TRATIO_CUT=0.5
SNR_CUT=10.0

# If this is can be corrected for in selection effect corrections, it is marked
# as such

MCAL_SELECT_FLAGS = {
    'dec_cut':2**0,
    'flags_badregion':2**1,
    'flags_gold':2**2,
    'mcal_flags':2**3,
    'size_cut':2**4,     # used in corrections
    'snr_cut':2**5,      # used in corrections
    'mof_flags':2**6,
    'mof_lowflux':2**7,
    'mcal_lowflux':2**8, # used in corrections
    'not_measured':FLAG_NOT_MEASURED,
}

# for these we also have the sheared versions for each
# name
MCAL_SHEAR_NAME_MAP = OrderedDict([
    ('e1',         {'name':'mcal_g',     'dt':'f8', 'el':0}),
    ('e2',         {'name':'mcal_g',     'dt':'f8', 'el':1}),
    ('size',       {'name':'mcal_T_r',   'dt':'f8', 'el':None}),
    ('size_err',   {'name':'mcal_T_err', 'dt':'f8', 'el':None}),
    ('snr',        {'name':'mcal_s2n_r', 'dt':'f8', 'el':None}),
    ('covmat_0_0', {'name':'mcal_g_cov', 'dt':'f8', 'el':(0,0)}),
    ('covmat_0_1', {'name':'mcal_g_cov', 'dt':'f8', 'el':(0,1)}),
    ('covmat_1_1', {'name':'mcal_g_cov', 'dt':'f8', 'el':(1,1)}),
])

MCAL_R_DT = [
    ('R11', 'f8'),
    ('R12', 'f8'),
    ('R21', 'f8'),
    ('R22', 'f8'),
]

MCAL_SHEAR_TYPES= [
    'noshear',
    '1p',
    '1m',
    '2p',
    '2m',
]

MOF_COLS_TO_COPY = [
    'ra',
    'dec',
    'tilename',
]
def get_nbands_from_pars(pars):
    return len(pars)-5

def get_bands(nbands):
    if nbands==3:
        bands = ['r','i','z']
    elif nbands==4:
        bands = ['g','r','i','z']
    else:
        raise ValueError("can't determine bands for nband==%d" % nbands)

    return bands

def get_mcal_name_map(nbands):

    bands=get_bands(nbands)

    nmap = copy.deepcopy(MCAL_NAME_MAP)

    for i,band in enumerate(bands):
        n=Namer(back=band)
        nmap[n('nimage_tot')] = {
            'name':'nimage_tot',
            'dt':'i4',
            'el':i,
        }
        nmap[n('nimage_use')] = {
            'name':'nimage_use',
            'dt':'i4',
            'el':i,
        }

    return nmap


def get_mcal_shear_name_map(nbands):

    bands=get_bands(nbands)

    nmap = copy.deepcopy(MCAL_SHEAR_NAME_MAP)
    for i,band in enumerate(bands):
        n=Namer(back=band)
        nmap[n('flux')] = {
            'name':'mcal_pars',
            'dt':'f8',
            'el':5+i,
        }
        nmap[n('flux_err')] = {
            'name':'mcal_pars_cov',
            'dt':'f8',
            'el':(5+i,5+i),
            'mcal_types':['noshear'],
        }

    return nmap

def set_defaults(st):
    """
    default for flag columns is always 2**31-1 "not_measured"
    """
    for dt in st.dtype.descr:
        name = dt[0]

        if 'cov' in name or 'err' in name:
            # default to large positive values for errors
            # and covariances
            st[name] = PDEFAULT

        elif 'S' in dt[1]:
            st[name] = SDEFAULT

        elif 'flag' in name.lower():
            st[name] = FLAG_NOT_MEASURED

        else:
            st[name] = DEFAULT


def load_mcal_flist(dir):
    pattern = os.path.join(dir, MCAL_FILE_PATTERN)

    flist = glob.glob(pattern)
    flist.sort()

    if len(flist) == 0:
        raise RuntimeError("found no files '%s'" % pattern)

    return flist

def get_namer(type=None):
    if type == 'noshear':
        back=None
    else:
        back=type

    return Namer(back=back)

class Namer(object):
    """
    create strings with a specified front prefix
    or back suffix
    """
    def __init__(self, front=None, back=None):
        if front=='':
            front=None
        if back=='':
            back=None

        self.front=front
        self.back=back

        if self.front is None and self.back is None:
            self.nomod=True
        else:
            self.nomod=False

    def __call__(self, name):
        if self.nomod:
            return name
        else:
            n=name
            if self.front is not None:
                n = '%s_%s' % (self.front, n)
            if self.back is not None:
                n = '%s_%s' % (n, self.back)
            return n



