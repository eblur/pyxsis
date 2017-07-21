import numpy as np
import clarsach

KEV  = ['kev', 'keV']
ANGS = ['Angstroms','Angstrom','Angs','angstroms','angstrom','angs']
ALLOWED_UNITS = KEV + ANGS

__all__ = ['Spectrum']

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

    def _parse_binned_counts(self):
        binning = self.binning[self.notice]
        return 0.0

    def _parse_binned_edges(self):
        binning = self.binning[self.notice]
        return 0.0

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
            counts = self._parse_binned_counts()
            ener_lo, ener_hi = self._parse_binned_edges()

        # Figure out how counts should be arranged
        if unit in KEV:
            sl = slice(None, None, 1)
            new_lo, new_hi = ener_lo, ener_hi
        if unit in ANGS:
            sl = slice(None, None, -1)
            new_lo = clarsach.respond.CONST_HC/ener_hi[sl]
            new_hi = clarsach.respond.CONST_HC/ener_lo[sl]

        new_mid = 0.5 * (new_lo + new_hi)
        return new_lo, new_hi, new_mid, noticed[sl]
