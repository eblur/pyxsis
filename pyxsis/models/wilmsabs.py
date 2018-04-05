import os
from astropy.table import Table
from scipy.interpolate import interp1d
from .model import Model

def _find_tablefile(name):
    root_path = os.path.dirname(__file__)
    data_path = root_path + '/tables/'
    return data_path + name

class WilmsAbs(Model):
    def __init__(self, nH=1.e22, lims=[(0, 1.e24)], units=['cm^-2']):
        Model.__init__(self, ['nH'], [nH], lims, units)
        self.xsect = self._read_xsect()
        return

    def _read_xsect(self):
        fname = _find_tablefile('cosplot_table.txt')
        tt = Table.read(fname, format='ascii')
        return interp1d(tt['col1'], tt['col2'])

    #def calculate(self, ener_lo, ener_hi):
