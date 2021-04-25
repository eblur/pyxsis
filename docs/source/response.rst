
Instrument Response Classes
===========================

Interpreting X-ray spectra requires knowledge of the instrumental
response. There are two essential types of response files. Pyxsis
provides a class for handling each of them.

Auxiliary Response File (ARF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Auxiliary Response File (frequenty referred to simply as **ARF**)
defines the telescope efficiency as a function of energy. This is
typically a 1D file with units of [cm^2 counts / photon] as a function
of photon energy. See the `CXC documentation page on ARFs
<https://cxc.cfa.harvard.edu/ciao/dictionary/arf.html>`_

.. autoclass:: pyxsis.xrayspectrum1d.ARF
    :members:

Redistribution Matrix File (RMF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Redistribution Matrix File (frequenty referred to simply as
**RMF**) describes the probability distribution of detector pulse
heights that arise when a photon interacts with the detector. It is a
2D matrix with the pulse height distribution as a function of incoming
photon energy. Seee the `CXC documentation page on RMFs <https://cxc.cfa.harvard.edu/ciao/dictionary/rmf.html>`_

.. autoclass:: pyxsis.xrayspectrum1d.RMF
    :members:
