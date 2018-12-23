## Plotting functions

import astropy.units as u
import numpy as np
from . import binspectrum

__all__ = ['plot_counts', 'plot_unfold', 'plot_model_flux']

def plot_counts(ax, spectrum, xunit='keV', perbin=True, rate=False, \
                bkgsub=True, usebackscal=True, **kwargs):

    lo, hi, cts, cts_err = spectrum.binned_counts(bkgsub=bkgsub, usebackscal=usebackscal)

    xlo = lo.to(u.Unit(xunit), equivalencies=u.spectral())
    xhi = hi.to(u.Unit(xunit), equivalencies=u.spectral())
    mid = 0.5 * (xlo + xhi)

    if rate:
        exp = spectrum.exposure
    else:
        exp = 1.0

    if perbin:
        dbin   = 1.0
    else:
        dbin   = np.abs(xhi - xlo)

    y    = cts/exp/dbin
    yerr = cts_err/exp/dbin

    ax.errorbar(mid.value, y.value, yerr=yerr.value,
                ls='', markersize=0, color='k', capsize=0, alpha=0.5)
    ax.step(mid, y, where='mid', **kwargs)
    ax.set_xlabel(mid.unit)
    ax.set_ylabel(y.unit)
    return


def plot_unfold(ax, spectrum, xunit='keV', perbin=False, \
                bkgsub=True, usebackscal=True, **kwargs):

    assert isinstance(spectrum, binspectrum.Spectrum)

    # Models will always be in keV bin units
    no_mod  = np.ones_like(spectrum.arf.specresp)  # a non-model of ones (integrated)
    eff_tmp = spectrum.apply_resp(no_mod)

    def _bin_exp(exp, binning):
        # Use mean effective exposure for binned spectra
        nstart, nend = min(binning), max(binning)
        result = [np.mean(exp[binning == n]) for n in np.arange(nstart, nend+1)]
        assert len(result) == (nend - nstart + 1)
        return np.array(result)

    # Now deal with desired xunit
    assert xunit in ALLOWED_UNITS
    lo, hi, cts, cts_err = spectrum.bin_counts(xunit, bkgsub=bkgsub, usebackscal=usebackscal)
    mid = 0.5 * (lo + hi)
    if all(spectrum.binning == 0):
        eff_exp = eff_tmp[spectrum.notice]
    else:
        eff_exp = _bin_exp(eff_tmp[spectrum.notice], spectrum.binning[spectrum.notice])

    flux, f_err = np.zeros_like(eff_exp), np.zeros_like(eff_exp)
    ii = np.isfinite(eff_exp) & (eff_exp != 0.0)
    if xunit in ANGS:
        flux[ii] = cts[ii] / eff_exp[ii][::-1]
        f_err[ii] = cts_err[ii] / eff_exp[ii][::-1]
    else:
        flux[ii] = cts[ii] / eff_exp[ii]
        f_err[ii] = cts_err[ii] / eff_exp[ii]

    if perbin:
        dbin   = 1.0
        ylabel = 'Flux [phot cm$^{-2}$ s$^{-1}$ bin$^{-1}$]'
    else:
        dbin   = hi - lo
        ylabel = 'Flux [phot cm$^{-2}$ s$^{-1}$ %s$^{-1}$]' % xunit

    # Now plot it
    ax.errorbar(mid, flux/dbin, yerr=f_err/dbin,
                ls='', marker=None, color='k', capsize=0, alpha=0.5)
    ax.step(lo, flux/dbin, where='mid', **kwargs)
    ax.set_xlabel("%s" % xunit)
    ax.set_ylabel(ylabel)

def plot_model_flux(ax, spectrum, model, xunit='keV', perbin=False, **kwargs):

    lo, hi, cts, cts_err = spectrum.bin_counts(unit=xunit)
    mid = 0.5 * (lo + hi)

    mflux = model.calculate(spectrum.arf.e_low, spectrum.arf.e_high)  # returns flux per bin
    if xunit in ANGS:
        mflux = mflux[::-1]

    if perbin:
        dbin   = 1.0
        ylabel = 'Flux [phot cm$^{-2}$ s$^{-1}$ bin$^{-1}$]'
    else:
        dbin   = hi - lo
        ylabel = 'Flux [phot cm$^{-2}$ s$^{-1}$ %s$^{-1}$]' % xunit

    ax.plot(mid, mflux/dbin, **kwargs)
    ax.set_xlabel("%s" % xunit)
    ax.set_ylabel(ylabel)
