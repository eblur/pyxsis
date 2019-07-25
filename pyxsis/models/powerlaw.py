from .model import Model
import numpy as np

__all__ = ['PowerLaw']

class PowerLaw(Model):
    def __init__(self, norm=1.e-10, phoindex=2.0,
                 lims=[(0, 1.e10), (-10.,10.)],
                 units=['s^-1 cm^-2', ''],
                 name='PowerLaw'):
        """
        PowerLaw Model
        --------------
        norm : float
            Normalization value (typically phot/s/cm^2)
        phoindex : float
            Photon index for the power law spectrum
        """
        Model.__init__(self, ['norm','phoindex'],
                       [norm, phoindex],
                       lims, units, name=name)
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
        Flux per bin (photon/cm**2/s)
        """
        assert len(ener_lo) == len(ener_hi)
        # Computes flux spectrum [phot cm^-2 s^-1] for given energy grid
        # integral over the power law model
        norm = self['norm']
        phoindex = self['phoindex']
        if phoindex == 1.0:
            r = norm * (np.log(ener_hi) - np.log(ener_lo))
        else:
            r = norm / (1.0-phoindex) * \
                ( ener_hi**(-phoindex + 1.0) - \
                  ener_lo**(-phoindex + 1.0) )
        return r
