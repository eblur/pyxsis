import numpy as np
import astropy.units as u
from astropy.io import fits
from .binspectrum import XBinSpectrum

''' NOT YET TESTED '''

KEV  = ['kev', 'keV']
ANGS = ['Angstroms','Angstrom','Angs','angstroms','angstrom','angs']
ALLOWED_UNITS = KEV + ANGS

ALLOWED_FORMATS = ['chandra_hetg','chandra_acis']

__all__ = ['XBkgSpectrum']

class XBkgSpectrum(XBinSpectrum):
    """
    Class for reading in background spectra. This is a subclass of XBinSpectrum.

    Parameters
    ----------
    bin_lo, bin_hi, counts, exposure, **kwargs (see also `specutils.XraySpectrum1D`)

    **To load from a file:**

    XBkgSpectrum(from_file=, format=, colname=)

    from_file : str
        Name of FITS file for background spectrum

    format : str
        Format of FITS file (Default: 'chandra_hetg' uses the 'BACKGROUND_UP'
        and 'BACKGROUND_DOWN' columns to read in the data)

    colname : str
        Name of FITS data column that contains the counts histogram of interest.
        (Default: 'COUNTS')

    Attributes
    ----------
    Inherits all attributes from XBinSpectrum

    **Additional attributes**

    backscale : numpy.ndarry or float
        Value for scaling the background count rate to the associated source area. 
        Defaults to 1.0
    """
    def __init__(self, *args, from_file=None, format='chandra_hetg', colname='COUNTS', **kwargs):
        if from_file is None:
            super().__init__(*args, **kwargs)
            self.backscale = 1.0
        else:
            if format == 'chandra_hetg':
                bin_lo, bin_hi, cts, exp, bscal = self._read_HETG(from_file)
                super().__init__(bin_lo, bin_hi, cts, exp, **kwargs)
                self.backscale = bscal
            else:
                bin_lo, bin_hi, cts, exp, bscal = self._read_other(from_file, colname)
                super().__init__(bin_lo, bin_hi, cts, exp, **kwargs)
                self.backscale = bscal

    def _read_HETG(self, filename):
        ff     = fits.open(filename)
        data   = ff[1].data
        hdr    = ff[1].header
        bin_unit = u.Unit(data.columns['BIN_LO'].unit)
        exposure = hdr['EXPOSURE'] * u.second
        counts = (data['BACKGROUND_UP'] + data['BACKGROUND_DOWN']) * u.ct
        bin_lo = data['BIN_LO'] * bin_unit
        bin_hi = data['BIN_HI'] * bin_unit
        # area of srouce region / area of background region
        backscal = hdr['BACKSCAL'] / (hdr['BACKSCUP'] + hdr['BACKSCDN'])
        ff.close()
        return bin_lo, bin_hi, counts, exposure, backscal

    def _read_other(self, filename, colname):
        ff     = fits.open(filename)
        data   = ff[1].data
        hdr    = ff[1].header
        bin_unit = u.Unit(data.columns['BIN_LO'].unit)
        exposure = hdr['EXPOSURE'] * u.second
        counts = data[colname] * u.ct
        bin_lo = data['BIN_LO'] * bin_unit
        bin_hi = data['BIN_HI'] * bin_unit
        # area of background region
        try:
            backscal = 1.0 / hdr['BACKSCAL']
        except:
            backscal = 1.0
        ff.close()
        return bin_lo, bin_hi, counts, exposure, backscal

    def bin_bkg(self, notice, binning, use_backscale=True):
        """
        Returns a binned background spectrum

        Parameters
        ----------
        notice : ndarray, dtype=bool
            Defines what regions of the spectrum to notice

        binning : ndarray
            Defines the binning for the spectrum (see Spectrum.binning)

        usebackscale : bool
            If True, the background will be scaled using XBkgSpectrum.backscale

        Returns
        -------
        bin_lo, bin_hi, bkg_counts, bkg_counts_err : ndarrays
        """
        # Deal with backscal, which could be an array
        backscal, scalar_backscal = 1.0, True
        if use_backscale:
            if np.size(self.backscale) == 1:
                backscal, scalar_backscal = self.backscale, True
            else:
                backscal, scalar_backscal = self.backscale[notice], False

        if all(binning == 0):
            bin_lo = self.bin_lo[notice]
            bin_hi = self.bin_hi[notice]
            result  = self.counts[notice] * backscal
            result_err = np.sqrt(self.counts[notice]) * backscal  # propogated error
        else:
            assert len(notice) == len(binning)  # need to apply notice array to binning
            binning = binning[notice]
            counts  = self.counts[notice]
            ener_lo = self.bin_lo[notice]
            ener_hi = self.bin_hi[notice]

            bin_lo, bin_hi, result, result_err2 = [], [], [], []
            for n in np.arange(min(binning), max(binning)+1):
                if scalar_backscal:
                    bb = backscal
                else:
                    bb = backscal[binning == n]
                bin_lo.append(ener_lo[binning == n][0])
                bin_hi.append(ener_hi[binning == n][-1])
                result.append(np.sum(counts[binning == n] * bb))
                result_err.append(np.sqrt(np.sum(counts[binning == n] * bb**2)))  # propogated error

            bin_lo = np.array(bin_lo)
            bin_hi = np.array(bin_hi)
            result = np.array(result)
            result_err = np.array(result_err)

        return bin_lo, bin_hi, result, result_err
