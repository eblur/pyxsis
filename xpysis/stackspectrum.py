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

        # Modify / overwrite old stuff
        self.counts       = self._stack_counts(speclist)
        self.exposure     = 1.0  # Pull from ARF, if exposure exists, if not throw a warning and take it from spectrum files
        self.arf.specresp = self._stack_arf(speclist)
        self.arf.exposure = 1.0
        self.arf.fracexpo = 1.0

    def _stack_counts(self, speclist):
        s0     = speclist[0]
        counts = np.copy(s0.counts)
        for s in speclist[1:]:
            assert all(s.bin_lo == s0.bin_lo), "Grids on every spectrum need to match"
            assert all(s.bin_hi == s0.bin_hi), "Grids on every spectrum need to match"
            assert s.bin_unit == s0.bin_unit, "Grids need to be in the same units"
            counts += s.counts
        return counts

    def _stack_arf(self, speclist):
        return 0.0
