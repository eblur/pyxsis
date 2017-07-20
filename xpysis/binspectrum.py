import numpy as np
import clarsach

KEV  = ['kev', 'keV']
ANGS = ['Angstroms','Angstrom','Angs','angstroms','angstrom','angs']
ALLOWED_UNITS = KEV + ANGS

class Spectrum(clarsach.XSpectrum):
    def __init__(self, filename, **kwargs):
        clarsach.XSpectrum.__init__(self, filename, **kwargs)
        self.hard_set_units('keV')  # Always keep binning in keV
        self.notice  = None
        self.binning = None

    @property
    def ener_mid(self):
        return 0.5 * (self.ener_lo + self.ener_hi)

    @property
    def angs_mid(self):
        return 0.5 * (self.angs_lo + self.angs_hi)

    def notice(self, bmin, bmax, unit='keV'):
        assert unit in ALLOWED_UNITS
        if unit in ANGS:
            emin = clarsach.CONST_HC / bmax
            emax = clarsach.CONST_HC / bmin
        if unit in KEV:
            emin, emax = bmin, bmax

        # Make sure the user hasn't changed bin_units to angstrom
        assert self.bin_unit in KEV
        self.notice = (self.bin_lo >= emin) & (self.bin_hi < emax)

    @property
    def bin_counts(self, unit='keV'):
        # Always make sure this is true
        assert self.bin_unit in KEV

        # Returns lo, hi, mid, counts, counts_err
        noticed = self.counts[self.notice]
        ener_lo = self.bin_lo[self.notice]
        ener_hi = self.bin_hi[self.notice]

        # Figure out how counts should be arranged
        if unit in KEV:
            sl = slice(None, None, 1)
            new_lo, new_hi = ener_lo, ener_hi
        if unit in ANGS:
            sl = slice(None, None, -1)
            new_lo = clarsach.CONST_HC/ener_hi[sl]
            new_hi = clarsach.CONST_HC/ener_lo[sl]

        new_mid = 0.5 * (new_lo, new_hi)
        return new_lo, new_hi, new_mid, noticed[sl], np.sqrt(noticed[sl])
