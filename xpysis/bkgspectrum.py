import os
import numpy as np
import clarsach

KEV  = ['kev', 'keV']
ANGS = ['Angstroms','Angstrom','Angs','angstroms','angstrom','angs']
ALLOWED_UNITS = KEV + ANGS

__all__ = ['BkgSpectrum']

class BkgSpectrum(object):
    """
    Class for reading in background spectra
    """
    def __init__(self, filename):
        self.__store_path(filename)

    def _read_chandra(self, filename):
        this_dir = os.path.dirname(os.path.abspath(filename))
        ff   = fits.open(filename)
        self.data = ff[1].data
        self.rmf_file = this_dir + "/" + ff[1].header['RESPFILE']
        self.arf_file = this_dir + "/" + ff[1].header['ANCRFILE']
        self.rmf = clarsach.RMF(self.rmf_file)
        self.arf = clarsach.ARF(self.arf_file)
        ff.close()
