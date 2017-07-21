import numpy as np
import clarsach

KEV  = ['kev', 'keV']
ANGS = ['Angstroms','Angstrom','Angs','angstroms','angstrom','angs']
ALLOWED_UNITS = KEV + ANGS

__all__ = ['Spectrum','group_channels']

class Spectrum(clarsach.XSpectrum):
    def __init__(self, filename, **kwargs):
        clarsach.XSpectrum.__init__(self, filename, **kwargs)
        self._setbins_to_keV()  # Always keep binning in keV
        self.notice  = np.ones_like(self.counts, dtype=bool)
        self.binning = np.zeros_like(self.counts)

    @property
    def ener_mid(self):
        return 0.5 * (self.ener_lo + self.ener_hi)

    @property
    def angs_mid(self):
        return 0.5 * (self.angs_lo + self.angs_hi)

    def notice_values(self, bmin, bmax, unit='keV'):
        assert unit in ALLOWED_UNITS
        if unit in ANGS:
            emin = clarsach.respond.CONST_HC / bmax
            emax = clarsach.respond.CONST_HC / bmin
        if unit in KEV:
            emin, emax = bmin, bmax

        # Make sure the user hasn't changed bin_units to angstrom
        assert self.bin_unit in KEV
        self.notice = (self.bin_lo >= emin) & (self.bin_hi < emax)

    def notice_all(self):
        # Resets the notice attribute
        self.notice = np.ones_like(self.counts, dtype=bool)

    def reset_binning(self):
        # Resets the binning attribute
        self.binning = np.zeros_like(self.counts)

    def _parse_binning(self):
        ## Returns the number of counts in each bin for a binned spectrum
        ## Works on the noticed regions only
        assert not all(self.binning == 0), "there is no grouping on this spectrum"

        # Use noticed regions only
        binning = self.binning[self.notice]
        counts  = self.counts[self.notice]
        ener_lo = self.bin_lo[self.notice]
        ener_hi = self.bin_hi[self.notice]

        bin_lo = np.array([ener_lo[binning == n][0] for n in np.arange(min(binning), max(binning)+1)])
        bin_hi = np.array([ener_hi[binning == n][-1] for n in np.arange(min(binning), max(binning)+1)])
        result = np.array([np.sum(counts[binning == n]) for n in np.arange(min(binning), max(binning)+1)])


        # Unit tests
        assert len(bin_lo) == (max(binning) - min(binning) + 1)
        assert len(bin_hi) == (max(binning) - min(binning) + 1)
        assert all(bin_lo < bin_hi)
        assert all(bin_lo[1:] == bin_hi[:-1])
        assert len(result) == (max(binning) - min(binning) + 1)
        assert np.sum(result) == np.sum(counts)  # Make sure no counts are lost

        return bin_lo, bin_hi, result

    def bin_counts(self, unit='keV'):
        # It's assumed that the spectrum is stored in keV bin units
        assert self.bin_unit in KEV

        # Returns lo, hi, mid, counts
        # self.binning is set to all zeros if there is no binning defined
        if all(self.binning == 0.0):
            counts  = self.counts[self.notice]
            ener_lo = self.bin_lo[self.notice]
            ener_hi = self.bin_hi[self.notice]
        else:
            ener_lo, ener_hi, counts = self._parse_binning()

        # Figure out how counts should be arranged
        if unit in KEV:
            sl = slice(None, None, 1)
            new_lo, new_hi = ener_lo, ener_hi
        if unit in ANGS:
            sl = slice(None, None, -1)
            new_lo = clarsach.respond.CONST_HC/ener_hi[sl]
            new_hi = clarsach.respond.CONST_HC/ener_lo[sl]

        new_mid = 0.5 * (new_lo + new_hi)
        return new_lo, new_hi, new_mid, counts[sl]

## ----- Binning functions

def group_channels(spectrum, n):
    """
    Group channels in a spectrum by a constant factor, n

    Parameters
    ----------
    spectrum : Spectrum
        Must contain `binning` attribute (ndarray)

    n : int
        Integer factor for binning the spectrum channels

    Returns
    -------
    Modifies Spectrum.binning, returns nothing
    """
    assert n > 1, "n must be larger than 1"

    tot_chan = len(spectrum.binning)  # Total number of channels
    counter, index = 0, 0
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
    spectrum : Spectrum
        Must contain `binning` attribute (ndarray)

    mc : int
        Minimum number of counts per bin

    Returns
    -------
    Modifies Spectrum.binning, returns nothing
    """
    assert mc > 0, "mc must be larger than 1"

    tot_chan = len(spectrum.binning)
    counter, tempcount = 0, 0
    result = []
    for i in range(tot_chan):
        result.append(counter)
        tempcount += spectrum.counts[i]
        if tempcount >= mc:
            counter += 1
            tempcount = 0

    spectrum.binning = np.array(result)
    return
