.. pyXsis documentation master file, created by
   sphinx-quickstart on Sat Mar 27 10:31:13 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the docmentation page for pyXsis
===========================================

pyXsis = Python X-ray Spectral Interpretation System

This is a toy python code for manipulating high resolution X-ray
Spectra.

This python library owes it's inspiration to the `Interactive Spectral
Interpretation System
<http://adsabs.harvard.edu/abs/2000ASPC..216..591H>`_ written
by J. Houck and J. Davis (creator of S-lang), as well as long-time
users and contributors M. Nowak, J. Wilms, and the group at Remeis
Observatory. Thank you also to the entire Chandra HETG group at MIT
for personal help throughout the years with interpreting high
resolution X-ray spectra.

This library also directly utilizes code from the `clarsach
<https://github.com/dhuppenkothen/clarsach>`_ repo written by Daniela
Huppenkothen, which reverse-engineered the response file treatment
provided by the standard high energy software package, `XSPEC
<https://heasarc.gsfc.nasa.gov/xanadu/xspec/>`_.

Installation
^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2
	      
   install

Package Basics
^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2
	      
   spectrum/index
   response/index
   io
   plotting
   examples/index


