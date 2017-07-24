import numpy as np
import clarsach
from binspectrum import Spectrum

KEV  = ['kev', 'keV']
ANGS = ['Angstroms','Angstrom','Angs','angstroms','angstrom','angs']
ALLOWED_UNITS = KEV + ANGS

__all__ = ['StackSpectrum']

class StackSpectrum(Spectrum):
    def __init__(self, file0, speclist, **kwargs):
        Spectrum.__init__(self, file0, **kwargs)

        assert isinstance(speclist, list), "Need to provide a list of Spectrum objects"
        assert len(speclist) > 1, "Need more than one spectrum to stack"
        self.counts   = self._stack_counts(self, speclist[0], speclist[1:])
        self.specresp = self._stack_arf(self, speclist[0], speclist[1:])
        self.exposure = None  # Pull from ARF, if exposure exists, if not throw a warning and take it from spectrum files
        self.fracexpo = 1.0
        
    def _stack_counts(self, s0, speclist):
        counts = s0.counts
        for s in speclist:
            assert all(s.bin_lo == s0.bin_lo), "Grids on every spectrum need to match"
            assert all(s.bin_hi == s0.bin_hi), "Grids on every spectrum need to match"
            assert s.bin_unit == s0.bin_unit, "Grids need to be in the same units"
            counts += s.counts
        return counts

    def _stack_arf(self, s0, speclist):
        return 0.0

    def apply_arf(self, mflux):
        return mflux * self.specresp

    def apply_resp(self, mflux, exposure=None):
        # Given a model flux spectrum, apply the response
        mrate = self.apply_arf(mflux)  # phot/s per bin
        if exposure is None:
            mrate *= self.exposure
        else:
            mrate *= exposure  # phot per bin

        result = self.rmf.apply_rmf(mrate)  # counts per bin
        return result

    @property
    def bin_mid(self):
        return 0.5 * (self.bin_lo + self.bin_hi)

    def _return_in_units(self, unit):
        assert unit in ALLOWED_UNITS
        if unit == self.bin_unit:
            return (self.bin_lo, self.bin_hi, self.bin_mid, self.counts)
        else:
            # Need to use reverse values if the bins are listed in increasing order
            new_lo, sl = clarsach.respond._Angs_keV(self.bin_hi)
            new_hi, sl = clarsach.respond._Angs_keV(self.bin_lo)
            new_mid = 0.5 * (new_lo + new_hi)
            new_cts = self.counts[sl]
            return (new_lo, new_hi, new_mid, new_cts)
