# pyxsis
_Python X-ray Spectral Interpretation System_

Toy python code for manipulating high resolution X-ray spectra

**Acknowledgements:**
This python library owes it's inspiration to the [Interactive Spectral Interpretation System](http://adsabs.harvard.edu/abs/2000ASPC..216..591H) written by J. Houck and J. Davis (creator of S-lang), as well as long-time users and contributors M. Nowak, J. Wilms, and the group at Remeis Observatory. Thank you also to the entire Chandra HETG group at MIT for personal help throughout the years with interpreting high resolution X-ray spectra.

## Install instructions

First download my version of _clarsach_ from Github and install it:

```
git clone https://github.com/eblur/clarsach.git
cd clarsach
python setup.py install
cd ..
```

Next, download the _pyxsis_ repository from Github and install it:

```
git clone https://github.com/eblur/pyxsis.git
cd pyxsis
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

## Installing a development version of pyXsis

These instructions teach you how to create an Anaconda environment for pyXsis development.

```
conda create -n pyxsis-dev python=3
source activate pyxsis-dev
conda install numpy scipy matplotlib astropy
```

Go to the folder where you would like to keep your libraries. Then install [Clarsach](https://github.com/dhuppenkothen/clarsach).
```
git clone git@github.com:eblur/clarsach.git
cd clarsach
python setup.py install
cd ..
```

Now clone and install pyXsis.
```
git clone git@github.com:eblur/pyxsis.git
cd pyxsis
python setup.py develop
```

If you use Jupyter notebooks, you can install the conda environment as a separate kernel, using the `ipykernel` package.
```
conda install ipykernel
python -m ipykernel install --user --name pyxsis-dev --display-name "pyxsis-dev"
```

When you are done playing with pyXsis, you can exit out the conda environment with:
```
source deactivate
```
