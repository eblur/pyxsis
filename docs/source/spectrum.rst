
1D Spectrum Classes
===================

This package provides several 1D spectrum classes and funcitons to
support common tasks used in visualizing and interpreting X-ray
spectra.

.. contents:: :local:
   :depth: 2


XraySpectrum1D
^^^^^^^^^^^^^^

The basic functionality to load and visualize X-ray spectra is
provided with the XraySpectrum1D class.

To initialize XraySpectrum1D, you must input the counts histogram data
directly.

.. autoclass:: pyxsis.xrayspectrum1d.XraySpectrum1D
    :members:

XBinSpectrum
^^^^^^^^^^^^

.. autoclass:: pyxsis.binspectrum.XBinSpectrum
    :members:

XBkgSpectrum
^^^^^^^^^^^^

.. autoclass:: pyxsis.binspectrum.XBkgSpectrum
    :members:
       
Binning Functions
^^^^^^^^^^^^^^^^^

Group by number of channels
---------------------------

.. autofunction:: pyxsis.binspectrum.group_channels

Group by number of counts
-------------------------

.. autofunction:: pyxsis.binspectrum.group_mincounts

		  
Apply a custom binning
----------------------

.. autofunction:: pyxsis.binspectrum.bin_anything
		  		 
