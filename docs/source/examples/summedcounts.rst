
Plot a summed counts histogram
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because it's nearly impossible to combine an RMF file from different
observations or instrument spectra (e.g., combining the RMF from a +1
and -1 order), **we do not recommend fitting combined spectra**. You
should use your preferred fitting method to fit separate spectra
simultaneously. However, in low signal-to-noise spectra it may be
useful to display the summed +1 and -1 orders from an HETG spectrum
for visual inspection only. This example demonstrates how to load two
spectra, sum their counts, create a new Spectrum object to hold the
summed counts histogram, bin the spectra, and plot the resulting count
rates.

**WARNING:** This method requires the two spectra to have the same bin
edge values, which is true when comparing +1 and -1 orders from the
same gratings type. No errors will arise if you attempt to use this
method to combine HEG+1 and MEG+1 orders, for example, but the results
will be entirely incorrect!

::

   import pyxsis

   # Load the HEG+1 and HEG-1 order spectra
   spec_hm1 = pyxsis.io.load_chandra_hetg("tests/data/17392/heg_-1.pha")
   spec_hp1 = pyxsis.io.load_chandra_hetg("tests/data/17392/heg_1.pha")

   # Sum the counts histograms
   summed_counts = spec_hm1.counts + spec_hp1.counts
   summed_heg = pyxsis.XBinSpectrum(spec_hm1.bin_lo, spec_hm1.bin_hi,
       summed_counts, spec_hm1.exposure)

   # Group the summed histogram
   # and apply that same grouping to the original spectra
   # (for visualization purposes)
   pyxsis.group_mincounts(summed_heg, 100)
   spec_hm1.binning = summed_heg.binning
   spec_hp1.binning = summed_heg.binning
   
   # Plot the count rates for each spectrum on the same axis
   ax = plt.subplot(111)
   ax.set_xlim(1.7,1.9)
   pyxsis.plot_counts(ax, spec_hm1, rate=True,
       color='r', label='HEG-1')
   pyxsis.plot_counts(ax, spec_hp1, rate=True,
       color='b', label='HEG+1')
   pyxsis.plot_counts(ax, summed_heg, rate=True,
       color='k', label='Summed')
   ax.legend()
   plt.show()


