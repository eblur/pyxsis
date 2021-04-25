Install Instructions
====================

Requirements
^^^^^^^^^^^^
* `Astropy <https://www.astropy.org/>`_
* `Specutils <https://specutils.readthedocs.io/en/stable/>`_
* Numpy
* Matplotlib

Install development version
^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can download and install pyXsis from the Github repository:
::

   git clone https://github.com/eblur/pyxsis.git
   cd pyxsis
   python setup.py install

Testing the Installation
========================

To run the test notebooks in pyxsis/tests/notebooks, you will need to
download the test HETG data.

First, download the tar ball from `<https://doi.org/10.5281/zenodo.2528474>`_

Next, copy the downloaded file to your pyxsis folder and unpack the tar ball:
::

   tar -xvf pyxsis_test_data.tar

You're now ready to run the notebooks!

Quick Start Guide
=================

To open a .pha file with pyXsis, you can follow these commands:
::

   import pyxsis
   my_spectrum = pyxsis.io.load_chandra_hetg('my_HETG_file.pha')

   import matplotlib.pyplot as plt
   ax = plt.subplot(111)
   pyxis.plot_counts(ax, my_spectrum, xunit='keV')

