## Plotting functions

import numpy as np
import clarsach

KEV      = ['kev', 'keV']
ANGS     = ['angs', 'angstrom', 'Angstrom', 'angstroms', 'Angstroms']

ALLOWED_UNITS = KEV + ANGS

__all__ = ['plot_counts', 'plot_unfold', 'plot_model_flux']

def plot_counts(ax, spectrum, xunit='keV', perbin=True, **kwargs):
    if isinstance(spectrum, clarsach.XSpectrum):
        lo, hi, mid, cts = spectrum._return_in_units(xunit)
    else:
        lo, hi, mid, cts = spectrum.bin_counts(xunit)
    cts_err = np.sqrt(cts)

    if perbin:
        dbin   = 1.0
        ylabel = 'Counts per bin'
    else:
        dbin   = hi - lo
        ylabel = 'Counts %s$^{-1}$' % xunit

    ax.errorbar(mid, cts/dbin, yerr=cts_err/dbin,
                ls='', markersize=0, color='k', capsize=0, alpha=0.5)
    ax.step(lo, cts/dbin, where='post', **kwargs)
    ax.set_xlabel("%s" % xunit)
    ax.set_ylabel(ylabel)
    return

def plot_unfold(ax, spectrum, xunit='keV', perbin=False, **kwargs):

    # Models will always be in keV bin units
    no_mod  = np.ones(len(spectrum.arf.specresp))  # a non-model of ones (integrated)
    eff_exp = spectrum.apply_resp(no_mod)  # non-counts per bin

    # Now deal with desired xunit
    assert xunit in ALLOWED_UNITS
    if isinstance(spectrum, clarsach.XSpectrum):
        lo, hi, mid, cts = spectrum._return_in_units(xunit)
    else:
        lo, hi, mid, cts = spectrum.bin_counts(xunit)

    flux, f_err = np.zeros_like(eff_exp), np.zeros_like(eff_exp)
    ii        = np.isfinite(eff_exp) & (eff_exp != 0.0)
    if xunit in ANGS:
        flux[ii] = cts[ii] / eff_exp[ii][::-1]
        f_err[ii] = np.sqrt(cts[ii]) / eff_exp[ii][::-1]
    else:
        flux[ii] = cts[ii] / eff_exp[ii]
        f_err[ii] = np.sqrt(cts[ii]) / eff_exp[ii]

    if perbin:
        dbin   = 1.0
        ylabel = 'Flux [phot cm$^{-2}$ s$^{-1}$ bin$^{-1}$]'
    else:
        dbin   = hi - lo
        ylabel = 'Flux [phot cm$^{-2}$ s$^{-1}$ %s$^{-1}$]' % xunit

    # Now plot it
    ax.errorbar(mid, flux/dbin, yerr=f_err/dbin,
                ls='', marker=None, color='k', capsize=0, alpha=0.5)
    ax.step(lo, flux/dbin, where='post', **kwargs)
    ax.set_xlabel("%s" % xunit)
    ax.set_ylabel(ylabel)

def plot_model_flux(ax, spectrum, model, xunit='keV', perbin=False, **kwargs):
    assert xunit in ALLOWED_UNITS

    if isinstance(spectrum, clarsach.XSpectrum):
        lo, hi, mid, cts = spectrum._return_in_units(xunit)
        elo, ehi, emid, cts = spectrum._return_in_units('keV')
    else:
        lo, hi, mid, cts = spectrum.bin_counts(xunit)
        elo, ehi, emid, cts = spectrum.bin_coutns('keV')

    mflux = model.calculate(spectrum.arf.ener_lo, spectrum.arf.ener_hi)  # returns flux per bin
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
