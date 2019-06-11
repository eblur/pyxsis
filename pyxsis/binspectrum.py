import numpy as np
import astropy.units as u

from .xrayspectrum1d import XraySpectrum1D
from .bkgspectrum import XBkgSpectrum

__all__ = ['XBinSpectrum','group_channels','group_mincounts','bin_anything']

class XBinSpectrum(XraySpectrum1D):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.notice  = np.ones_like(self.counts.value, dtype=bool)
        self.binning = np.zeros_like(self.counts.value, dtype=int)
        self.bkg = None
        self.model = None
        self.model_counts = None

    @staticmethod
    def load(filename, format='chandra_hetg', arf=None, rmf=None):
        temp   = XraySpectrum1D.read(filename, format=format)
        result = XBinSpectrum(temp.bin_lo, temp.bin_hi,
                              temp.counts, temp.exposure)
        if arf is not None:
            result.assign_arf(arf)

        if rmf is not None:
            result.assign_rmf(rmf)

        return result

    ##-- I wrote these properties for convenience
    @property
    def bmid_keV(self):
        return self.spectral_axis.to(u.keV, equivalencies=u.spectral())

    @property
    def bmid_angs(self):
        return self.spectral_axis.to(u.angstrom, equivalencies=u.spectral())

    def notice_range(self, bmin, bmax):
        """
        Define edges for spectral regions to notice. Notices regions exclusively.

        >>> XBinSpectrum.notice(1.0*u.keV, 3.0*u.keV)

        will notice *only* the region between 1 and 3 keV.

        Parameters
        ----------
        bmin : astropy.Quantity
            Minimum value for the notice region

        bmax : astropy.Quantity
            Maxium value for the notice region

        Returns
        -------
        Modifies the XraySpectrum1D.notice attribute.
        """
        bmin0 = bmin.to(self.bin_lo.unit, equivalencies=u.spectral())
        bmax0 = bmax.to(self.bin_lo.unit, equivalencies=u.spectral())
        imin  = (self.bin_lo >= min(bmin0, bmax0))
        imax  = (self.bin_hi <= max(bmin0, bmax0))
        self.notice  = (imin & imax)

    def notice_all(self):
        """
        Resets the notice attribute to notice the entire range.
        """
        self.notice = np.ones_like(self.counts, dtype=bool)

    def reset_binning(self):
        """
        Resets the binning attribute
        """
        self.binning = np.zeros_like(self.counts)

    def binned_counts(self, subtract_bkg=False, use_backscale=True):
        """
        Returns the binned counts histogram from the noticed spectral region.

        Parameters
        ----------
        subtract_bkg : bool
            If True, supply the background subtracted region spectrum
            (only if a background spectrum is supplied)

        use_backscale : bool
            If True, scales the background by the backscal value

        Returns
        -------
        bin_lo : astropy.Quantity
            Lower edges for the new bins

        bin_hi : astropy.Quantity
            Higher edges for the new bins

        counts : astropy.Quantity
            Counts for the new bins

        cts_err : astropy.Quantity
            Error on the new bins
        """
        if all(self.binning == 0.0):
            counts  = self.counts[self.notice]
            bin_lo = self.bin_lo[self.notice]
            bin_hi = self.bin_hi[self.notice]
            cts_err = np.sqrt(counts.value) * u.ct
        else:
            bin_lo, bin_hi, counts, cts_err = self._parse_binning()

        if subtract_bkg and (self.bkg is not None):
            blo, bhi, bcts, bcts_err = self.binned_bkg(use_backscale=use_backscale)
            new_counts = counts - bcts
            new_error  = np.sqrt(cts_err**2 + bcts_err**2)
            return blo, bhi, new_counts, new_error
        else:
            return bin_lo, bin_hi, counts, cts_err

    def _parse_binning(self):
        ## Returns the number of counts in each bin for a binned spectrum
        ## Works on the noticed regions only
        assert not all(self.binning == 0), "there is no grouping on this spectrum"

        # Use noticed regions only
        binning = self.binning[self.notice]
        counts  = self.counts[self.notice].value
        bin_lo  = self.bin_lo[self.notice].value
        bin_hi  = self.bin_hi[self.notice].value

        new_bin_lo = np.array([bin_lo[binning == n][0] for n in np.arange(min(binning), max(binning)+1)])
        new_bin_hi = np.array([bin_hi[binning == n][-1] for n in np.arange(min(binning), max(binning)+1)])
        result = np.array([np.sum(counts[binning == n]) for n in np.arange(min(binning), max(binning)+1)])

        # Make sure no counts are lost
        percent_diff = np.abs((np.sum(result) - np.sum(counts))/np.sum(counts))
        assert percent_diff < 1.e-5, print("Binned counts are off by {} \%".format(percent_diff/100.))

        new_bin_lo *= self.bin_lo.unit
        new_bin_hi *= self.bin_lo.unit
        result *= u.ct
        result_err = np.sqrt(result.value) * u.ct

        return new_bin_lo, new_bin_hi, result, result_err

    def binned_bkg(self, use_backscale=True):
        """
        Return the background spectrum using the same spectral binning.

        Parameters
        ----------
        use_backscale : bool
            If True, returns the background scaled by the backscal value.

        Returns
        -------
        bin_lo : astropy.Quantity
            Lower edges for the new background bins

        bin_hi : astropy.Quantity
            Higher edges for the new background bins

        counts : astropy.Quantity
            Counts for the binned background histogram

        cts_err : astropy.Quantity
            Error on the new background bins
        """
        return self.bkg.binned_counts(self.notice, self.binning, use_backscale=use_backscale)

    def assign_bkg(self, bkgspec, format='chandra_hetg', colname='COUNTS'):
        """
        Assign a background spectrum.

        Parameters
        ----------
        bkgspec : str or XBkgSpectrum

            If a string, loads background from the specified FITS file.
            Otherwise, stores the input in the XBinSpectrum.bkg attributes.

        format : str

            Specifies the format of the background file.

            If 'chandra_hetg', runs `XBkgSpectrum.load_HETG`. Otherwise, runs
            `XBkgSpectrum.load`

        colname : str

            Specifies the column name of the relevant counts histogram when
            running `XBkgSpectrum.load`
        """
        if isinstance(bkgspec, str):
            if format == 'chandra_hetg':
                self.bkg = XBkgSpectrum.load_HETG(bkgspec)
            else:
                self.bkg = XBkgSpectrum.load(bkgspec, colname=colname)
        else:
            self.bkg = bkgspec

    def evaluate_model(self, model_flux):
        model_w_arf = self.arf.apply_arf(model_flux)
        model_w_rmf = self.rmf.apply_rmf(model_w_arf)
        self.model_counts = model_w_rmf

    def bin_model_counts(self, unit='keV'):
        # Use noticed regions only
        binning = self.binning[self.notice]
        counts  = self.model_counts[self.notice]
        result  = np.array([np.sum(counts[binning == n]) for n in np.arange(min(binning), max(binning)+1)])
        # Figure out how counts should be arranged
        if unit in KEV:
            sl = slice(None, None, 1)
        if unit in ANGS:
            sl = slice(None, None, -1)
        return result[sl]

## ----- Binning functions

def group_channels(spectrum, n):
    """
    Group channels in a spectrum by a constant factor, n

    Parameters
    ----------
    spectrum : pyxsis.XBinSpectrum
        Must contain `binning` attribute (ndarray)

    n : int
        Integer factor for binning the spectrum channels

    Returns
    -------
    Modifies spectrum.binning
    """
    assert n > 1, "n must be larger than 1"

    tot_chan = len(spectrum.binning)  # Total number of channels
    counter, index = 1, 0
    result = np.array([])
    while index < tot_chan:
        result = np.append(result, np.repeat(counter, n))
        index += n
        counter += 1

    # Trim array in case we went over the total number of channels
    # which will happen when tot_chan % n != 0
    result = result[np.arange(tot_chan)]
    spectrum.binning = result
    return

def group_mincounts(spectrum, mc):
    """
    Group channels in a spectrum so that there is a minimum number of counts in each bin

    Parameters
    ----------
    spectrum : XBinSpectrum
        Must contain `binning` attribute (ndarray)

    mc : int
        Minimum number of counts per bin

    Returns
    -------
    Modifies spectrum.binning
    """
    assert mc > 0, "mc must be larger than 1"

    tot_chan = len(spectrum.binning)
    counter, tempcount = 0, 0
    result = []
    for i in range(tot_chan):
        result.append(counter)
        tempcount += spectrum.counts[i].value
        if tempcount >= mc:
            counter += 1
            tempcount = 0
    result = np.array(result)

    # Ensure that the last bin has at least ten counts
    nlast = result[-1]
    last_bin_count = np.sum(spectrum.counts[result == nlast].value)
    if last_bin_count < mc:
        result[result == nlast] = nlast-1

    spectrum.binning = result
    return

def bin_anything(x, binning, notice=None):
    """
    Group anything according to a binning array

    Parameters
    ----------
    x : np.array
        Array to bin

    binning : np.array
        Bin numbers for sorting

    notice : bool np.array (optional)
        Array values to notice

    Returns
    -------
    A numpy array that holds the binned results
    """
    assert len(x) == len(binning)
    if notice is None:
        notice = np.ones_like(x, dtype='bool')
    xx = x[notice]
    bb = binning[notice]
    if np.all(bb == 0):
        return xx
    else:
        return np.array([np.sum(xx[bb == n]) for n in np.arange(min(bb), max(bb)+1)])

