from .model import Model
import numpy as np

__all__ = ['PowerLaw']

class PowerLaw(Model):
    def __init__(self, norm=1.e-10, phoindex=2.0,
                 lims=[(0, 1.e10), (-10.,10.)],
                 units=['s^{-1} cm^{-2}', '']):
        """
        PowerLaw Model
        --------------
        norm : Normalization value (typically phot/s/cm^2)
        phoindex : Photon index ()
        """
        Model.__init__(self, ['norm','phoindex'], [norm, phoindex], lims, units)
        return
