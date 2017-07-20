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
    def bin_counts(self):
        noticed = self.counts[self.notice]
        return noticed
