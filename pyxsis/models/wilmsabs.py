import os
import numpy as np
import astropy.units as u
from astropy.table import Table
from scipy.interpolate import interp1d
from .model import Model

def _find_tablefile(name):
    root_path = os.path.dirname(__file__)
    data_path = root_path + '/tables/'
    return data_path + name

class WilmsAbs(Model):
    def __init__(self, nH=1.e22, lims=[(0, 1.e24)], units=['cm^-2'],
                 name='WilmsAbs'):
        """
        WilmsAbs Model
        --------------
        Wilms et al. (2000) model for interstellar absorption.
        Uses table from tbnew model, as of 2015.07.08

        nH : float
            Hydrogen column density [cm^-2]
        """
        Model.__init__(self, ['nH'], [nH], lims, units, name=name)
        self.xsect = self._read_xsect()
        return

    def _read_xsect(self):
        """
        Reads tabulated model from Wilms et al. (2000)

        Returns
        -------
        scipy.interpolate.interp1d object of energy (keV) vs extinction (cm^2)
        """
        fname = _find_tablefile('cosplot_table.txt')
        tt = Table.read(fname, format='ascii')
        return interp1d(tt['col1'], tt['col2'])

    def calculate(self, ener_lo, ener_hi):
        """
        Calculates the extinction factor, exp(-tau) for ISM absorption

        Parameters
        ----------
        ener_lo : numpy.ndarray
            Low energy edge of counts histogram

        ener_hi : numpy.ndarray
            High energy edge of counts histogram

        Returns
        -------
        exp(-tau_abs)
        """
        emid = 0.5 * (ener_lo + ener_hi)
        tau  = self.xsect(emid) * u.cm**2 * self['nH']
        return np.exp(-tau)
