from .model import Model
import numpy as np

__all__ = ['PowerLaw']

class PowerLaw(Model):
    def __init__(self, norm=1.e-10, phoindex=2.0,
                 lims=[(0, 1.e10), (-10.,10.)],
                 units=['s^-1 cm^-2', '']):
        """
        PowerLaw Model
        --------------
        norm : float
            Normalization value (typically phot/s/cm^2)
        phoindex : float
            Photon index for the power law spectrum
        """
        Model.__init__(self, ['norm','phoindex'], [norm, phoindex], lims, units)
        return

    ## Modified version of clarsach models
    def calculate(self, ener_lo, ener_hi):
        """
        Calculates the photon flux spectrum [phot cm^-2 s^-1]

        Parameters
        ----------
        ener_lo : numpy.ndarray
            Low energy edge of counts histogram

        ener_hi : numpy.ndarray
            High energy edge of counts histogram

        Returns
        -------
        Flux = norm * np.power(emid, -phoindex) where emid = 0.5 * (ener_lo + ener_hi) and has units of keV
        """
        assert len(ener_lo) == len(ener_hi)
        # Computes flux spectrum [phot cm^-2 s^-1] for given energy grid
        # integral over the power law model
        norm = self['norm'].value
        phoindex = self['phoindex'].value
        if phoindex == 1.0:
            r = np.log(ener_hi) - np.log(ener_lo)
        else:
            r = -norm * ener_hi**(-phoindex + 1.0) + \
                norm * ener_lo**(-phoindex + 1.0)
        return r
