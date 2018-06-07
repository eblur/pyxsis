from .model import Model
import numpy as np
import astropy.units as u

__all__ = ['Gaussian']

class Gaussian(Model):
    def __init__(self, norm=1.e-5, ecenter=6.7, sigma=0.1,
                 lims=[(0, 1.e10), (0.0, 200.0), (0.0,100.0)],
                 units=['s^-1 cm^-2', 'keV', 'keV'],
                 name='Gaussian'):
        """
        PowerLaw Model
        --------------
        norm : float
            Normalization value (typically phot/s/cm^2)

        ecenter : float (keV)
            Center value for the Gaussian

        sigma : float (keV)
            Width of the Gaussian
        """
        Model.__init__(self, ['norm','ecenter','sigma'],
                       [norm, ecenter, sigma],
                       lims, units, name=name)
        return

    ## Modified version of clarsach models
    def calculate(self, ener_lo, ener_hi):
        """
        Calculates the Gaussian spectrum in flux [phot cm^-2 s^-1]

        Parameters
        ----------
        ener_lo : numpy.ndarray
            Low energy edge of counts histogram

        ener_hi : numpy.ndarray
            High energy edge of counts histogram

        Returns
        -------
        Flux = norm * exp( -0.5 * (emid - ecenter)**2 / sigma**2 )
            where emid = 0.5 * (ener_lo + ener_hi) and has units of keV
        """
        assert len(ener_lo) == len(ener_hi)
        # Computes flux spectrum [phot cm^-2 s^-1] for given energy grid
        # integral over the power law model
        norm    = self['norm']
        ecenter = self['ecenter']
        sigma   = self['sigma']
        emid    = 0.5 * (ener_lo + ener_hi) * u.keV
        return norm * np.exp( -0.5 * (emid - ecenter)**2 / sigma**2 )
