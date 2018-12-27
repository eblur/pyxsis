import numpy as np
from specutils.xrayspectrum import XraySpectrum1D, AreaResponse, ResponseMatrix
from .binspectrum import XBinSpectrum

__all__ = ['stack_spectra']

def stack_spectra(speclist, rmf=None):
    """
    Stack a list of Spectrum objects

    Parameters
    ----------
    speclist : list of specutils.xrayspectrum.XraySpectrum1D objects
        A list of spectra to be stacked.

    rmf : specutils.xrayspectrum.ResponseMatrix object or string (file name)
        Response matrix to assign to the output spectrum.
        If _None_, the RMF from the first spectrum in speclist will be assigned
        to the final result.

    Returns
    -------
    XBinSpectrum with the stacked counts histogram, attached to a modified ARF
    and an RMF.
    """
    assert len(speclist) > 1, "Need more than one spectrum to stack"

    def _stack_counts(speclist):
        s0     = speclist[0]
        counts = np.copy(s0.counts)
        for s in speclist[1:]:
            assert all(s.bin_lo == s0.bin_lo), "Grids on every spectrum need to match"
            assert all(s.bin_hi == s0.bin_hi), "Grids on every spectrum need to match"
            counts += s.counts
        return counts

    def _arf_exp(spec):
        if spec.arf.exposure is None:
            print("Cannot find exposure keyword from ARF, using PHA file")
            return spec.exposure
        else:
            return spec.arf.exposure

    def _stack_arf(speclist):
        # Do a time weighted average of the ARF response and fracexpo
        exposure = _arf_exp(speclist[0])
        specresp = speclist[0].arf.eff_area.value * exposure
        fracexpo = speclist[0].arf.fracexpo * exposure
        elow0    = speclist[0].arf.e_low
        ehigh0   = speclist[0].arf.e_high

        for s in speclist[1:]:
            a, exp = s.arf, _arf_exp(s)
            assert all(a.e_low == elow0), "Grids on every arf need to match"
            assert all(a.e_high == ehigh0), "Grids on every arf need to match"
            specresp += a.specresp * exp
            exposure += exp
            fracexpo += a.fracexpo

        specresp /= exposure
        fracexpo /= exposure

        return elow0, ehigh0, specresp, fracexpo, exposure

    # Modify / overwrite old stuff
    new_counts = _stack_counts(speclist)
    elo, ehigh, arfresp, fracexpo, exposure = _stack_arf(speclist)
    new_arf = AreaResponse(elo, ehigh, arfres, fracexpo=fracexpo, exposure=exposure)
    return XBinSpectrum(speclist[0].bin_lo, speclist[0].bin_hi, new_counts, exposure,
                           arf=new_arf, rmf=speclist[0].rmf)
