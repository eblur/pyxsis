import numpy as np
from specutils import AreaResponse, ResponseMatrix
from .binspectrum import XBinSpectrum

__all__ = ['stack_spectra']

def stack_spectra(speclist, rmf=None, sum_exposure=False):
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

    sum_exposure : bool
        If _True_, the exposure times from each spectrum will be summed and the
        effective area will be computed with a time weighted average. Use this
        setting when combining spectra from different observations.
        If _False_, the exposure time will not be summed, and the effective area
        will be computed by summing the ARFs.

    Returns
    -------
    XBinSpectrum with the stacked counts histogram, with a modified ARF
    and an RMF.
    """
    assert len(speclist) > 1, "Need more than one spectrum to stack"

    # method for stacking counts
    def _stack_counts(speclist):
        s0     = speclist[0]
        counts = np.copy(s0.counts) * s0.counts.unit
        for s in speclist[1:]:
            assert all(s.bin_lo == s0.bin_lo), "Grids on every spectrum need to match"
            assert all(s.bin_hi == s0.bin_hi), "Grids on every spectrum need to match"
            counts += s.counts
        return counts

    # method for grabbing the exposure time
    def _arf_exp(spec):
        if spec.arf.exposure is None:
            print("Cannot find exposure keyword from ARF, using PHA file")
            return np.copy(spec.exposure) * spec.exposure.unit
        else:
            return np.copy(spec.arf.exposure) * spec.arf.exposure.unit

    # method for calculating time averaged ARF
    def _time_avg_arf(speclist):
        exposure = _arf_exp(speclist[0])
        arf0     = speclist[0].arf
        specresp = np.copy(arf0.eff_area) * arf0.eff_area.unit * exposure
        fracexpo = np.copy(arf0.fracexpo) * exposure
        elow0    = np.copy(arf0.e_low) * arf0.e_low.unit
        ehigh0   = np.copy(arf0.e_high) * arf0.e_high.unit

        for s in speclist[1:]:
            a, exp = s.arf, _arf_exp(s)
            assert all(a.e_low == elow0), "Grids on every arf need to match"
            assert all(a.e_high == ehigh0), "Grids on every arf need to match"
            specresp += a.eff_area * exp
            exposure += exp
            fracexpo += a.fracexpo * exp

        specresp /= exposure
        fracexpo /= exposure

        return elow0, ehigh0, specresp, fracexpo, exposure

    # method for calculating a summed arf
    def _summed_arf(speclist):
        exposure = _arf_exp(speclist[0])
        arf0     = speclist[0].arf
        specresp = np.copy(arf0.eff_area) * arf0.eff_area.unit
        fracexpo = np.copy(arf0.fracexpo)
        elow0    = np.copy(arf0.e_low) * arf0.e_low.unit
        ehigh0   = np.copy(arf0.e_high) * arf0.e_high.unit

        for s in speclist[1:]:
            a, exp = s.arf, _arf_exp(s)
            assert all(a.e_low == elow0), "Grids on every arf need to match"
            assert all(a.e_high == ehigh0), "Grids on every arf need to match"
            assert (exp == exposure), \
                "These spectra don't appear to be from the same observation." \
                "Run stack_spectra with sum_exposure = True."
            specresp += a.eff_area
            fracexpo += a.fracexpo

        fracexpo /= len(speclist)
        return elow0, ehigh0, specresp, fracexpo, exposure

    # Run all the internal functions
    new_counts = _stack_counts(speclist)
    if sum_exposure:
        elo, ehigh, specresp, fracexpo, exposure = _time_avg_arf(speclist)
    else:
        elo, ehigh, specresp, fracexpo, exposure = _summed_arf(speclist)

    # Create a new arf
    new_arf = AreaResponse(elo, ehigh, specresp, fracexpo=fracexpo, exposure=exposure)

    if rmf is None:
        rmf = speclist[0].rmf

    # Return the new XBinSpectrum with new arf and rmf from first spectrum attached
    return XBinSpectrum(speclist[0].bin_lo, speclist[0].bin_hi, new_counts, exposure,
                           arf=new_arf, rmf=rmf)
