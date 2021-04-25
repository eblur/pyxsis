
Plotting Functions
==================

This package provides convenience plotting functions for visualizing
and interpreting X-ray spectra. It relies on the `matplotlib
<https://matplotlib.org/>`_ plotting infrastructure.

All functions take a `Matplotlib AxesSubplot object
<https://matplotlib.org/stable/api/axes_api.html#the-axes-class>`_
object as the primary argument. We recommend setting up the Axes
object using the `matplotlib.pyplot.subplot
<https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplot.html>`_
function. For example, to set up axes for a single plot window:
::

   import matplotlib.pyplot as plt
   ax = plt.subplot(111)


Plot a counts histogram
^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: pyxsis.plot.plot_counts


Plot a flux spectrum
^^^^^^^^^^^^^^^^^^^^

.. autofunction:: pyxsis.plot.plot_unfold

