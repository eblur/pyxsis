import os
import numpy as np
import astropy.units as u
from .model import Model

## rgscat.py -- Provide a Rayleigh-Gans + Drude approximation scattering
## model for high imaging resolution spectral observations (e.g. Chandra and XMM)

class RGscat(Model):
    def __init__(self, tau1=0.5, lims=[(0., 100.)], units=[''],
                 name='RGscat'):
        """
        Rayleigh-Gans Scattering model
        ----------------------------------------------------------
        Uses the Rayleigh-Gans + Drude approximation to ISM dust scattering
        in the X-ray regime (see Mauche & Gorenstein 1986). Assumes that
        all of the dust scattering halo is removed from the spectral
        extraction aperture, and the scattering cross-section follows E^-2.

        tau1 : float
            Optical depth to scattering at 1 keV [unitless]
        """
        Model.__init__(self, ['tau1'], [tau1], lims, units, name=name)
        return

    def calculate(self, ener_lo, ener_hi):
        """
        Calculates the extinction factor, exp(-tau) for ISM dust scattering.
        This is a very rough calculation! Not to be used on highly binned data.

        Parameters
        ----------
        ener_lo : numpy.ndarray
            Low energy edge of counts histogram [keV]

        ener_hi : numpy.ndarray
            High energy edge of counts histogram [keV]

        Returns
        -------
        exp(-tau_sca(E)) where E = 0.5 * (ener_lo + ener_hi)
        """
        emid = 0.5 * (ener_lo + ener_hi)
        tau  = self['tau1'] / emid**2
        return np.exp(-tau)
