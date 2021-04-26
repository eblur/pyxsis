# pyXsis
_Python X-ray Spectral Interpretation System_

Toy python code for manipulating high resolution X-ray spectra. The documentation webpages are available here: https://pyxsis.readthedocs.io/en/latest/index.html 

**Acknowledgements:**

This python library owes it's inspiration to the [Interactive Spectral Interpretation System](http://adsabs.harvard.edu/abs/2000ASPC..216..591H) written by J. Houck and J. Davis (creator of S-lang), as well as long-time users and contributors M. Nowak, J. Wilms, and the group at Remeis Observatory. Thank you also to the entire Chandra HETG group at MIT for personal help throughout the years with interpreting high resolution X-ray spectra.

This library also directly utilizes code from the [clarsach](https://github.com/dhuppenkothen/clarsach) repo written by Daniela Huppenkothen, which reverse-engineered the response file treatment provided by the standard high energy software package, [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/).

## Install instructions


You must have [Astropy](http://www.astropy.org/) and [Specutils](https://specutils.readthedocs.io/en/latest/) installed.

Download and install _pyXsis_ from the Github repository:
```
git clone https://github.com/eblur/pyxsis.git
cd pyxsis
python setup.py install
```

### Testing the installation

To run the test notebooks in _pyxsis/tests/notebooks_, you will need to download the test HETG data.

First, download the tar ball from https://doi.org/10.5281/zenodo.2528474

Next, copy the downloaded file to your _pyxsis_ folder and unpack the tar ball:
```
tar -xvf pyxsis_test_data.tar
```

You're now ready to run the notebooks!

## Dependencies

+ Numpy Version 1.19
+ Astropy Version 4.2
+ Specutils Version 1.2
+ Matplotlib Version 3.3.4

## Quick start

```
import pyxsis
my_spectrum = pyxsis.io.load_chandra_hetg('my_Chandra_HETG_file.pha')

import matplotlib.pyplot as plt
ax = plt.subplot(111)
pyxsis.plot_counts(ax, my_spectrum, xunit='keV')
```

## Installing a development version of pyXsis

These instructions teach you how to create an Anaconda environment for pyXsis development.

```
conda create -n pyxsis-dev python=3
source activate pyxsis-dev
conda install numpy scipy matplotlib astropy specutils
```

Now clone and install pyXsis.
```
git clone https://github.com/eblur/pyxsis.git
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
