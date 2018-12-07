import numpy as np
import astropy.units as u
from specutils import XraySpectrum1D
#from .bkgspectrum import BkgSpectrum

KEV  = ['kev', 'keV']
ANGS = ['Angstroms','Angstrom','Angs','angstroms','angstrom','angs']
ALLOWED_UNITS = KEV + ANGS

#__all__ = ['Spectrum','group_channels','group_mincounts']

class XBinSpectrum(XraySpectrum1D):
    def __init__(self, from_file=None, format='chandra_hetg', **kwargs):
        if from_file is None:
            # this should never really happen. Need to figure out work around for testing
            XraySpectrum1D.__init__(self, np.array([]), np.array([]), u.angstrom,
                                    []*u.ct, 1.0*u.second, **kwargs)
        else:
            XraySpectrum1D.read(self, from_file, format=format, **kwargs)
        self.notice  = np.ones_like(self.counts, dtype=bool)
        self.binning = np.zeros_like(self.counts)
        self.bkg = None

    # Need to fix this for the fact that the lo and hi edges differ depending on unit
    def notice_values(self, bmin, bmax, unit='keV'):
        bin_edges  = np.append(self.bin_lo, self.bin_hi[-1])
        unit_edges = self.bin_edges.to(u.Unit(unit), equivalencies=u.spectral())
        self.notice = (unit_edges >= bmin) & (unit_edges <= bmax)[:-1]



        assert unit in ALLOWED_UNITS
        if unit in ANGS:
            emin = clarsach.respond.CONST_HC / bmax
            emax = clarsach.respond.CONST_HC / bmin
        if unit in KEV:
            emin, emax = bmin, bmax

        # Make sure the user hasn't changed bin_units to angstrom
        assert self.bin_unit in KEV
        self.notice = (self.bin_lo >= emin) & (self.bin_hi < emax)

'''
    def notice_all(self):
        # Resets the notice attribute
        self.notice = np.ones_like(self.counts, dtype=bool)

    def reset_binning(self):
        # Resets the binning attribute
        self.binning = np.zeros_like(self.counts)

    def bin_counts(self, unit='keV', bkgsub=True, usebackscal=True):
        # It's assumed that the spectrum is stored in keV bin units
        assert self.bin_unit in KEV

        # Returns lo, hi, mid, counts
        # self.binning is set to all zeros if there is no binning defined
        if all(self.binning == 0.0):
            counts  = self.counts[self.notice]
            ener_lo = self.bin_lo[self.notice]
            ener_hi = self.bin_hi[self.notice]
            cts_err = np.sqrt(counts)
        else:
            ener_lo, ener_hi, counts, cts_err = self._parse_binning()

        # Figure out how counts should be arranged
        if unit in KEV:
            sl = slice(None, None, 1)
            new_lo, new_hi = ener_lo, ener_hi
        if unit in ANGS:
            sl = slice(None, None, -1)
            new_lo = clarsach.respond.CONST_HC/ener_hi[sl]
            new_hi = clarsach.respond.CONST_HC/ener_lo[sl]

        if bkgsub and (self.bkg is not None):
            blo, bhi, bcts, bcts_err = self.bkg.bin_bkg(self.notice, self.binning, usebackscal=usebackscal)
            new_counts = counts[sl] - bcts[sl]
            new_error  = np.sqrt(cts_err[sl]**2 + bcts_err[sl]**2)
            return new_lo, new_hi, new_counts, new_error
        else:
            return new_lo, new_hi, counts[sl], np.sqrt(counts[sl])

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

        return bin_lo, bin_hi, result, np.sqrt(result)

    def bin_bkg(self, usebackscal=True):
        # Shortcut function for returning background
        return self.bkg.bin_bkg(self.notice, self.binning, usebackscal=usebackscal)

    def assign_bkg(self, bkgspec):
        assert isinstance(bkgspec, BkgSpectrum)
        assert all(self.bin_lo == bkgspec.bin_lo), "Grids need to be the same"
        assert all(self.bin_hi == bkgspec.bin_hi), "Grids need to be the same"
        assert self.bin_unit == bkgspec.bin_unit, "Bin units need to be the same"
        self.bkg = bkgspec

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
    result = np.array(result)

    # Ensure that the last bin has at least ten counts
    nlast = result[-1]
    last_bin_count = np.sum(spectrum.counts[result == nlast])
    if last_bin_count < mc:
        result[result == nlast] = nlast-1

    spectrum.binning = result
    return
'''
