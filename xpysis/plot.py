## Plotting functions

import numpy as np

ALLOWED_UNITS      = ['keV','angs','angstrom','kev']
CONST_HC    = 12.398418573430595   # Copied from ISIS, [keV angs]
UNIT_LABELS = dict(zip(ALLOWED_UNITS, ['Energy (keV)', 'Wavelength (angs)']))

__all__ = ['plot_counts','plot_unfold','plot_model_flux']

def _get_counts(spectrum, xunit='keV', **kwargs):
    lo, hi, mid, cts = spectrum._change_units(xunit)
    cts_err = np.sqrt(cts)
    # Return plot data if the user wants it
    return lo, hi, mid, cts, cts_err

def plot_counts(ax, spectrum, xunit='keV', perbin=True, **kwargs):
    lo, hi, mid, cts, cts_err = _get_counts(spectrum, xunit=xunit, **kwargs)

    if perbin:
        dbin   = 1.0
        ylabel = 'Counts per bin'
    else:
        dbin   = hi - lo
        ylabel = 'Counts %s$^{-1}$' % xunit

    ax.errorbar(mid, cts/dbin, yerr=cts_err/dbin,
                ls='', marker=None, color='k', capsize=0, alpha=0.5)
    ax.step(lo, cts/dbin, where='post', **kwargs)
    ax.set_xlabel(UNIT_LABELS[xunit])
    ax.set_ylabel(ylabel)

def plot_unfold(ax, spectrum, xunit='keV', perbin=False, **kwargs):
    eff_exp  = spectrum._eff_exposure()  # cm^2 sec count phot^-1
    # Have to take account of zero values in effective exposure
    flux, f_err = np.zeros(len(eff_exp)), np.zeros(len(eff_exp))
    ii        = (eff_exp != 0.0)
    flux[ii]  = spectrum.counts[ii] / eff_exp[ii]
    f_err[ii] = np.sqrt(spectrum.counts[ii]) / eff_exp[ii]

    # Now deal with desired xunit
    lo, hi, mid, cts = spectrum._change_units(xunit)
    if spectrum.bin_unit != xunit:
        flx, fe = flux[::-1], f_err[::-1]
    else:
        flx, fe = flux, f_err

    if perbin:
        dbin   = 1.0
        ylabel = 'Flux [phot cm$^{-2}$ s$^{-1}$ bin$^{-1}$]'
    else:
        dbin   = hi - lo
        ylabel = 'Flux [phot cm$^{-2}$ s$^{-1}$ %s$^{-1}$]' % xunit

    # Now plot it
    ax.errorbar(mid, flx/dbin, yerr=fe/dbin,
                ls='', marker=None, color='k', capsize=0, alpha=0.5)
    ax.step(lo, flx/dbin, where='post', **kwargs)
    ax.set_xlabel(UNIT_LABELS[xunit])
    ax.set_ylabel(ylabel)

def plot_model_flux(ax, spectrum, model, xunit='keV', perbin=False, **kwargs):
    lo, hi, mid, cts = spectrum._change_units(xunit)
    elo, ehi, emid, cts = spectrum._change_units('keV')
    mflux = model.calculate(elo, ehi)
    if xunit in ['angs','Angs','angstrom','Angstrom']:
        mflux = mflux[::-1]

    if perbin:
        dbin   = 1.0
        ylabel = 'Flux [phot cm$^{-2}$ s$^{-1}$ bin$^{-1}$]'
    else:
        dbin   = hi - lo
        ylabel = 'Flux [phot cm$^{-2}$ s$^{-1}$ %s$^{-1}$]' % xunit

    ax.plot(mid, mflux/dbin, **kwargs)
    ax.set_xlabel(UNIT_LABELS[xunit])
    ax.set_ylabel(ylabel)
