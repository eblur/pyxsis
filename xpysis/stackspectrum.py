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

        def _stack_counts(speclist):
            s0     = speclist[0]
            counts = np.copy(s0.counts)
            for s in speclist[1:]:
                assert all(s.bin_lo == s0.bin_lo), "Grids on every spectrum need to match"
                assert all(s.bin_hi == s0.bin_hi), "Grids on every spectrum need to match"
                assert s.bin_unit == s0.bin_unit, "Grids need to be in the same units"
                counts += s.counts
            return counts

        def _arf_exp(spec):
            result = 1.0
            try:
                result = spec.arf.data[1].header['EXPOSURE']
            except:
                result = spec.exposure
                print("Cannot find exposure keyword from ARF, using PHA file")
            return result

        def _stack_arf(speclist):
            # Do a time weighted average of the ARF response and fracexpo
            a0       = speclist[0].arf
            exposure = _arf_exp(speclist[0])
            specresp = np.copy(a0.specresp) * exposure
            fracexpo = np.copy(a0.fracexpo) * exposure

            for s in speclist[1:]:
                a, exp = s.arf, _arf_exp(s)
                assert all(a.e_low == a0.e_low), "Grids on every arf need to match"
                assert all(a.e_high == a0.e_high), "Grids on every arf need to match"
                assert a.e_unit == a0.e_unit, "ARF grids need to be in the same units"
                specresp += a.specresp * exp
                exposure += exp
                fracexpo += a.fracexpo

            specresp /= exposure
            fracexpo /= exposure

            return specresp, fracexpo, exposure
        
        # Modify / overwrite old stuff
        self.counts       = _stack_counts(speclist)
        arfresp, fracexpo, exposure = _stack_arf(speclist)
        self.arf.specresp = arfresp
        self.arf.fracexpo = fracexpo
        self.exposure     = exposure  # Taken from ARF if available
