## Plotting functions

import astropy.units as u
import numpy as np
from . import binspectrum

__all__ = ['plot_counts', 'plot_unfold']

def plot_counts(ax, spectrum, xunit='keV', perbin=True, rate=False,
                plot_bkg=False, subtract_bkg=True, use_backscale=True,
                scale_factor=1.0, **kwargs):
    """
    Plots the counts histogram for a 1D X-ray spectrum.

    ax : matplotlib AxesSubplot object
        The figure axis on which to plot.

    spectrum : pyxsis XBinSpectrum object
        The binnable 1D spectrum object to plot.

    xunit : string (default: 'keV')
        Defines the unit type to be used on the x-axis. The options
        are 'angs' and 'keV'.

    perbin : bool (default: True)
        If True, the plot y-axis will show the number of counts per
        bin. If False, the plot y-axis will show the number of counts
        per bin-width, depending on the xunit. For example, if
        xunit='keV', then the y-axis units will be counts/keV.

    rate : bool (default: False)
        If True, the plot y-axis will show the number of counts per
        second. If False, it will show the total number of counts.

    plot_bkg : bool (default: False)
        If True, the background spectrum assigned to the input
        spectrum will be plotted instead of the primary source
        spectrum.

    subtract_bkg : bool (default: True)
        If True, the background spectrum assigned to the input
        spectrum will be subtracted before plotting. If False, no
        background subtraction will be implemented.

    use_backscale : bool (default: True)
        If True, the background spectrum will be scaled by the
        pyxsis.XBkgSpectrum.backscale attribute. This attribute
        generally holds the ratio of the background extraction area to
        the source extraction area. If False, the raw background
        spectrum count rate will be used. This is helpful if you just
        want to view the raw background spectrum.

    scale_factor : float (default: 1.0)
        A normalization value to apply to the entire spectrum to be
        plotted. This option is provided solely for convenience, for
        example, in comparing spectra to each other.

    *kwargs* are passed to the main histrogram plotting function
    (ax.step). This can be used to change the color, line widths, line
    style, and more.
    """
    if plot_bkg:
        lo, hi, cts, cts_err = spectrum.binned_bkg(bin_unit=xunit,
                                                   use_backscale=use_backscale)
    else:
        lo, hi, cts, cts_err = spectrum.binned_counts(bin_unit=xunit,
                                                      subtract_bkg=subtract_bkg,
                                                      use_backscale=use_backscale)
    mid = 0.5 * (lo + hi)

    if rate:
        exp = spectrum.exposure
    else:
        exp = 1.0

    if perbin:
        dbin   = 1.0
    else:
        dbin   = np.abs(hi - lo)

    y    = cts/exp/dbin * scale_factor
    yerr = cts_err/exp/dbin * scale_factor

    ax.errorbar(mid.value, y.value, yerr=yerr.value,
                ls='', markersize=0, color='k', capsize=0, alpha=0.5)
    ax.step(mid, y, where='mid', **kwargs)
    ax.set_xlabel(mid.unit.to_string(format='latex_inline'))
    ax.set_ylabel(y.unit.to_string(format='latex_inline'))
    return


def plot_unfold(ax, spectrum, xunit='keV', perbin=False,
                subtract_bkg=True, use_backscale=True,
                scale_factor=1.0, **kwargs):
    """
    Plots the flux histogram for a 1D X-ray spectrum.

    ax : matplotlib AxesSubplot object
        The figure axis on which to plot.

    spectrum : pyxsis XBinSpectrum object
        The binnable 1D spectrum object to plot.

    xunit : string (default: 'keV')
        Defines the unit type to be used on the x-axis. The options
        are 'angs' and 'keV'.

    perbin : bool (default: False)
        If True, the plot y-axis will show the number of counts per
        bin. If False, the plot y-axis will show the number of counts
        per bin-width, depending on the xunit. For example, if
        xunit='keV', then the y-axis units will be counts/keV.

    subtract_bkg : bool (default: True)
        If True, the background spectrum assigned to the input
        spectrum will be subtracted before plotting. If False, no
        background subtraction will be implemented.

    use_backscale : bool (default: True)
        If True, the background spectrum will be scaled by the
        pyxsis.XBkgSpectrum.backscale attribute. This attribute
        generally holds the ratio of the background extraction area to
        the source extraction area. If False, the raw background
        spectrum count rate will be used. This is helpful if you just
        want to view the raw background spectrum.

    scale_factor : float (default: 1.0)
        A normalization value to apply to the entire spectrum to be
        plotted. This option is provided solely for convenience, for
        example, in comparing spectra to each other.

    *kwargs* are passed to the main histrogram plotting function
    (ax.step). This can be used to change the color, line widths, line
    style, and more.
    """
    # Models will always be in keV bin units
    # a non-model of ones (integrated)
    # This set of calcs will get the effective response for each bin (cm^2 ct / phot)
    no_mod  = np.ones(len(spectrum.counts))
    eff_tmp = spectrum.apply_response(no_mod)

    def _bin_exp(exp, binning):
        # Use mean effective exposure for binned spectra
        nstart, nend = min(binning), max(binning)
        result = [np.mean(exp[binning == n]) for n in np.arange(nstart, nend+1)]
        assert len(result) == (nend - nstart + 1)
        return np.array(result)

    # Get the binned counts
    lo, hi, cts, cts_err = spectrum.binned_counts(bin_unit=xunit,
                                                  subtract_bkg=subtract_bkg,
                                                  use_backscale=use_backscale)

    mid = 0.5 * (lo + hi)

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
        dbin   = np.abs(hi - lo)[ii]

    # plot values
    x = mid[ii]
    y = flux / dbin * scale_factor
    yerr = ferr / dbin * scale_factor

    # Now plot it
    ax.errorbar(x.value, y.value, yerr=yerr.value,
                ls='', marker=None, color='k', capsize=0, alpha=0.5)
    ax.step(x, y, where='mid', **kwargs)
    ax.set_xlabel(x.unit.to_string(format='latex_inline'))
    ax.set_ylabel("phot {}".format(y.unit.to_string(format='latex_inline')))

## Not yet tested
'''def plot_model_flux(ax, spectrum, model, scale_factor=1.0,
                       xunit='keV', perbin=False, **kwargs):

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

    ax.plot(mid, mflux/dbin * scale_factor, **kwargs)
    ax.set_xlabel("%s" % xunit)
    ax.set_ylabel(ylabel)'''
