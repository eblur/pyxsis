
Telescope Specific Loaders
==========================

The `pyxsis.io` submodule provides helper functions for loading files
that are specific to an X-ray observatory.


Chandra X-ray Observatory
^^^^^^^^^^^^^^^^^^^^^^^^^

To load a Chandra PHA level 1 file extracted from an HETG observation:
::
   
   from pyxsis.io import load_chandra_hetg
   spectrum = load_chandra_hetg("my_heg_m1.pha")

If the PHA file has the ARF and RMF filenames specified in the FITS
file header ('ANCRFILE' and 'RESPFILE', respectively), then they will
be automatically loaded by the pyxsis ARF and RMF classes. Otherwise,
no response files will be assigned and the user must specify them.
   
.. autofunction:: pyxsis.io.load_chandra_hetg
