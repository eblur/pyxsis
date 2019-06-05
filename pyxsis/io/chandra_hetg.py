from astropy.io import fits
from astropy.units import Unit

from .. import XraySpectrum1D

__all__ = ['load_chandra_hetg']

def load_chandra_hetg(file_name, arf=None, rmf=None):
    """
    Load Chandra HETG spectral data from a file into a spectrum object.

    Parameters
    ----------
    file_name: str
        The path to the FITS file

    arf : str OR ARF
        Filename for the area response file (ARF) or a pre-loaded AreaResponse object

    rmf : str OR RMF
        Filename for the response matrix file (RMF) or a pre-loaded ResponseMatrix object

    Returns
    -------
    data: XraySpectrum1D
        The spectrum that is represented by the data in this table.
    """

    with fits.open(file_name) as hdu:
        header = hdu[0].header
        meta   = {'header': header}
        data   = hdu[1].data

        bin_unit = Unit(data.columns['BIN_LO'].unit)
        bin_lo   = data['BIN_LO'] * bin_unit
        bin_hi   = data['BIN_HI'] * bin_unit

        counts   = data['COUNTS'] * Unit('ct')
        exposure = hdu[1].header['EXPOSURE'] * Unit('second')

    return XraySpectrum1D(bin_lo, bin_hi, counts,
                          exposure=exposure, arf=arf, rmf=rmf)
