# pyxsis
_Python X-ray Spectral Interpretation System_

Toy python code for manipulating high resolution X-ray spectra

**Acknowledgements:**
This python library owes it's inspiration to the [Interactive Spectral Interpretation System](http://adsabs.harvard.edu/abs/2000ASPC..216..591H) written by J. Houck and J. Davis (creator of S-lang), as well as long-time users and contributors M. Nowak, J. Wilms, and the group at Remeis Observatory. Thank you also to the entire Chandra HETG group at MIT for personal help throughout the years with interpreting high resolution X-ray spectra.

## Install instructions

Download the repository from github

```
git clone https://github.com/eblur/xpysis.git
```

Move into the _xpysis_ directory and run _setup.py_

```
python setup.py install
```

## Dependencies

+ Numpy Version 1.1 or later
+ Astropy Version 2.0 or later
+ [Clarsach](https://github.com/dhuppenkothen/clarsach) Version 0.0 or later

## Quick start

```
import pyxsis
my_spectrum = pyxsis.Spectrum('my_Chandra_HETG_file.pha', telescope='HETG')

import matplotlib.pyplot as plt
ax = plt.subplot(111)
pyxsis.plot_counts(ax, my_spectrum, xunit='kev')
```
