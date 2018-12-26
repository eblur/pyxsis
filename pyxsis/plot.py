## Plotting functions

import astropy.units as u
import numpy as np
from . import binspectrum

__all__ = ['plot_counts', 'plot_unfold']

def plot_counts(ax, spectrum, xunit='keV', perbin=True, rate=False, \
                plot_bkg=False, subtract_bkg=True, use_backscale=True, **kwargs):

    if plot_bkg:
        lo, hi, cts, cts_err = spectrum.binned_bkg(use_backscale=use_backscale)
    else:
        lo, hi, cts, cts_err = spectrum.binned_counts(subtract_bkg=subtract_bkg, use_backscale=use_backscale)

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
    ax.set_xlabel(mid.unit.to_string(format='latex_inline'))
    ax.set_ylabel(y.unit.to_string(format='latex_inline'))
    return


def plot_unfold(ax, spectrum, xunit='keV', perbin=False, \
                subtract_bkg=True, use_backscale=True, **kwargs):

    # Models will always be in keV bin units
    no_mod  = np.ones_like(spectrum.arf.specresp)  # a non-model of ones (integrated)
    eff_tmp = spectrum.apply_response(no_mod)

    def _bin_exp(exp, binning):
        # Use mean effective exposure for binned spectra
        nstart, nend = min(binning), max(binning)
        result = [np.mean(exp[binning == n]) for n in np.arange(nstart, nend+1)]
        assert len(result) == (nend - nstart + 1)
        return np.array(result)

    # Get the binned counts
    lo, hi, cts, cts_err = spectrum.binned_counts(subtract_bkg=subtract_bkg, use_backscale=use_backscale)

    xlo = lo.to(u.Unit(xunit), equivalencies=u.spectral())
    xhi = hi.to(u.Unit(xunit), equivalencies=u.spectral())
    mid = 0.5 * (xlo + xhi)

    # Get the binned effective area
    if all(spectrum.binning == 0):
        eff_exp = eff_tmp[spectrum.notice]
    else:
        eff_exp = _bin_exp(eff_tmp[spectrum.notice], spectrum.binning[spectrum.notice])
    eff_exp *= u.cm**2 * u.ct # cm^2 ct / phot

    # Calculate photon flux from binned effective area
    #flux, f_err = np.zeros_like(eff_exp), np.zeros_like(eff_exp)
    ii = np.isfinite(eff_exp) & (eff_exp != 0.0)
    flux = cts[ii] / eff_exp[ii] / spectrum.exposure # phot / cm^2 / second
    ferr = np.sqrt(cts[ii].value)*u.ct / eff_exp[ii] / spectrum.exposure

    if perbin:
        dbin   = 1.0
    else:
        dbin   = np.abs(xhi - xlo)[ii]

    # plot values
    x = mid[ii]
    y = flux / dbin
    yerr = ferr / dbin

    # Now plot it
    ax.errorbar(x.value, y.value, yerr=yerr.value,
                ls='', marker=None, color='k', capsize=0, alpha=0.5)
    ax.step(x, y, where='mid', **kwargs)
    ax.set_xlabel(x.unit.to_string(format='latex_inline'))
    ax.set_ylabel("phot {}".format(y.unit.to_string(format='latex_inline')))

## Not yet tested
'''def plot_model_flux(ax, spectrum, model, xunit='keV', perbin=False, **kwargs):

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
    ax.set_ylabel(ylabel)'''
