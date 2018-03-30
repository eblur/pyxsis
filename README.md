# xpysis
Python code for manipulating high resolution X-ray spectra

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
import xpysis
my_spectrum = xpysis.Spectrum('my_Chandra_HETG_file.pha', telescope='HETG')

import matplotlib.pyplot as plt
ax = plt.subplot(111)
xpysis.plot_counts(ax, my_spectrum, xunit='kev')
```
