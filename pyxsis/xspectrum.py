import os
import numpy as np
from astropy.io import fits
import astropy.units as u
#from clarsach.respond import RMF, ARF, _Angs_keV
from specutils import Spectrum1D

__all__ = ['XSpectrum']

KEV      = ['kev','keV']
ANGS     = ['angs','angstrom','Angstrom','angstroms','Angstroms']

ALLOWED_UNITS = KEV + ANGS
ALLOWED_TELESCOPES = ['HETG','ACIS']

def _read_chandra(filename):
    """
    Inputs
    ------
    pha file filename

    Returns
    -------
    bin_lo
    bin_hi
    bin_unit
    counts
    rmf_file
    arf_file
    """
    this_dir = os.path.dirname(os.path.abspath(filename))
    ff   = fits.open(filename)
    data = ff[1].data
    bin_lo   = data['BIN_LO']
    bin_hi   = data['BIN_HI']
    bin_unit = data.columns['BIN_LO'].unit
    counts   = data['COUNTS']
    rmf_file = this_dir + "/" + ff[1].header['RESPFILE']
    arf_file = this_dir + "/" + ff[1].header['ANCRFILE']
    exposure = ff[1].header['EXPOSURE']  # seconds
    ff.close()
    return bin_lo, bin_hi, bin_unit, counts, rmf_file, arf_file, exposure


# This code is from my dev version of clarsach
# https://github.com/eblur/clarsach/blob/master/clarsach/spectrum.py
class XSpectrum(Spectrum1D):
    def __init__(self, filename, telescope='HETG'):
        """
        Inputs
        ------
        filename : string
            Name of pha file to open

        telescope : string ['HETG' | 'ACIS']
            String description of instrument to use,
            only Chandra is supported right now

        Attributes
        ----------
        bin_lo : numpy.array
            Low end of the bin edges

        bin_hi : numpy.array
            High end of the bin edges

        bin_unit : string
            Description of the bin unit (see ALLOWED_UNITS)

        exposure : float
            Exposure time for the pha file (seconds)

        bin_unit : string ['angs' | 'kev']
            Description of the bin unit (see ALLOWED_UNITS)

        arf : clarsarch.respond.ARF object
            If an arf was not specified or found, will be None

        rmf : clarsach.respond.RMF object
            If an rmf was not specified or found, will be None

        """
        assert telescope in ALLOWED_TELESCOPES
        # Right now only Chandra is supported
        if telescope == 'HETG':
            bin_lo, bin_hi, bin_unit, counts, rmf, arf, exposure = _read_chandra(filename)
        elif telescope == 'ACIS':
            bin_lo, bin_hi, bin_unit, counts, rmf, arf, exposure = _read_chandra(filename)

        # instantiate with Spectrum1D
        bin_mid = 0.5 * (bin_lo + bin_hi)
        spectral_unit = ''
        if bin_unit in ANGS:
            spectral_unit = u.angstrom
        if bin_unit in KEV:
            spectral_unit = u.keV
        print("Spectral unit is {}".format(spectral_unit))
        Spectrum1D.__init__(self, spectral_axis=bin_mid, flux=counts,
                            spectral_axis_unit=spectral_unit, unit='')

        # Add other attributes
        self.bin_lo   = bin_lo
        self.bin_hi   = bin_hi
        self.bin_unit = bin_unit
        self.exposure = exposure
        print("loading rmf file: {}".format(rmf))
        self.rmf = RMF(rmf)
        print("loading arf file: {}".format(arf))
        self.arf = ARF(arf)
        self.__store_path(filename)

        if self.bin_unit != self.arf.e_unit:
            print("Warning: ARF units and pha file units are not the same!!!")
        if self.bin_unit != self.rmf.energ_unit:
            print("Warning: RMF units and pha file units are not the same!!!")
        return

    @property
    def counts(self):
        return self.flux

    def __store_path(self, filename):
        self.path = '/'.join(filename.split('/')[0:-1]) + "/"
        return

    def apply_resp(self, mflux, exposure=None):
        """
        Given a model flux spectrum, apply the response. In cases where the
        spectrum has both an ARF and an RMF, apply both. Otherwise, apply
        whatever response is in RMF.
        The model flux spectrum *must* be created using the same units and
        bins as in the ARF (where the ARF exists)!
        Parameters
        ----------
        mflux : iterable
            A list or array with the model flux values in ergs/keV/s/cm^-2
        exposure : float, default None
            By default, the exposure stored in the ARF will be used to compute
            the total counts per bin over the effective observation time.
            In cases where this might be incorrect (e.g. for simulated spectra
            where the pha file might have a different exposure value than the
            ARF), this keyword provides the functionality to override the
            default behaviour and manually set the exposure time to use.
        Returns
        -------
        count_model : numpy.ndarray
            The model spectrum in units of counts/bin
        """

        if self.arf is not None:
            mrate  = self.arf.apply_arf(mflux, exposure=exposure)
        else:
            mrate = mflux

        count_model = self.rmf.apply_rmf(mrate)

        return count_model

    @property
    def bin_mid(self):
        return 0.5 * (self.bin_lo + self.bin_hi)

    def _return_in_units(self, unit):
        assert unit in ALLOWED_UNITS
        if unit == self.bin_unit:
            return (self.bin_lo, self.bin_hi, self.bin_mid, self.counts)
        else:
            # Need to use reverse values if the bins are listed in increasing order
            new_lo, sl = _Angs_keV(self.bin_hi)
            new_hi, sl = _Angs_keV(self.bin_lo)
            new_mid = 0.5 * (new_lo + new_hi)
            new_cts = self.counts[sl]
            return (new_lo, new_hi, new_mid, new_cts)

    def _setbins_to_keV(self):
        assert self.bin_unit in ANGS
        new_bhi, sl = _Angs_keV(self.bin_lo)
        new_blo, sl = _Angs_keV(self.bin_hi)
        new_cts  = self.counts[sl]

        # Now hard set everything
        self.bin_lo = new_blo
        self.bin_hi = new_bhi
        self.flux   = new_cts
        self.bin_unit = si.keV
        return
