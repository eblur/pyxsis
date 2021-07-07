
Model Classes
=============

The pyXsis `Model` superclass provides basic functions for storing
parameter names, values, units, and limits.

Every the parameter value can be accessed by calling the Model object
like a dictionary. For example, to see the default photon index value
for the PowerLaw model:
::
	
	from pyxsis.models import PowerLaw
	my_pl = PowerLaw()
	my_pl['phoindex']


.. toctree::

   examples


Class Docstrings
^^^^^^^^^^^^^^^^

.. autoclass:: pyxsis.models.Model
    :members:


