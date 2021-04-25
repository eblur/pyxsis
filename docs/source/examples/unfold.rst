
Group a spectrum and plot it
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This example shows how to load a spectrum with pyXsis, group the
spectrum to a minimum number of counts, and plot the resulting flux
spectrum.
::

   import matplotlib.pyplot as plt
   from pyxsis.io import load_chandra_hetg
   from pyxsis import group_mincounts, plot_unfold

   spec = load_chandra_hetg("tests/data/17385/heg_-1.pha")
   group_mincounts(spec, 30)

   ax = plt.subplot(111)
   plot_unfold(ax, spec, xunit='keV')
   plt.loglog()
   plt.show()

