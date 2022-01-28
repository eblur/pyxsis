import os
import numpy as np

from astropy.io import fits
from astropy.units import Unit

from .. import XBinSpectrum

__all__ = ['load_xmm_rgs']

def load_xmm_rgs(filename, arf=None, rmf=None):
    """
    Load XMM RGS spectral data from a file into a spectrum object.

    **Inputs**
    
    filename : str
        The path to the FITS file

    arf : str -or- pyxsis.xrayspectrum1d.ARF
        Filename for the area response file (ARF) or a pre-loaded AreaResponse object
        For RGS data, there seems to be no separate ARF file. So this will usually be None

    rmf : str -or- pyxsis.xrayspectrum1d.RMF
        Filename for the response matrix file (RMF) or a pre-loaded ResponseMatrix object

    **Returns**
    
    pyxsis XraySpectrum1D object representing the data in the input FITS file
    """
    this_dir = os.path.dirname(os.path.abspath(filename))

    with fits.open(filename) as hdu:
        header  = hdu[0].header
        meta    = {'header': header}
        data    = hdu[1].data
        datahdr = hdu[1].header

        bin_unit  = Unit(datahdr['TCUNI1'])
        bin_cen   = datahdr['TCRVL1'] # center of the reference bin
        bin_cen_i = datahdr['TLMIN1']-1 # reference bin (assuming 1 refers to first bin)
        bin_width = datahdr['TCDLT1'] # width of each channel
        bin_lo    = (np.arange(datahdr['TLMAX1']) * bin_width + bin_cen - bin_width / 2.0) * bin_unit
        bin_hi    = (np.arange(datahdr['TLMAX1']) * bin_width + bin_cen + bin_width / 2.0) * bin_unit

        counts   = data['COUNTS'] * Unit('count')
        exposure = datahdr['EXPOSURE'] * Unit('second')

        if rmf is None:
            rmf = this_dir + "/" + hdu[1].header['RESPFILE']

    return XBinSpectrum(bin_lo, bin_hi, counts,
                          exposure=exposure, arf=arf, rmf=rmf)
