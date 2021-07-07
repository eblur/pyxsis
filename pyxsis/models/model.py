import astropy.units as u

class Model(object):
    """Superclass for handling X-ray spectral models
    
    **Inputs**
    
    par_keys : list of strings : Parameter names
    
    par_vals : list of floats : Initial parameter values
    
    par_lims : list of tuples : Lists of parameter limits
    
    par_units : list of strings : Unit string for each parameter
    
    name : string : A name describing the model
    
    **Attributes**
    
    name
    
    keys
    
    vals : dictionary of key-value pairs for model parameters
    
    lims : dictionary of key-tuple pairs for model paraemters
    
    units : dictionary of key-unit pairs for model parameters
    """
    def __init__(self, par_keys, par_vals, par_lims, par_units, name='Model'):
        # Check that the parameter values fall within the limits
        # If they don't, return None
        for v,l in zip(par_vals, par_lims):
            assert self._check_lims(v, l)
            assert l[1] > l[0]
        self.name = name
        self.keys = par_keys
        self.vals = dict(zip(par_keys, par_vals))
        self.lims = dict(zip(par_keys, par_lims))
        self.units = dict(zip(par_keys, par_units))
        return

    def _check_lims(self, val, lim):
        """
        **Inputs**
        
        val : float
        lim : tuple (floats)

        **Returns**

        True if `val` is within the boundaries of `lim`
        """
        result = (val >= lim[0]) & (val <= lim[1])
        return result

    def check_par_lims(self, new_par_dict):
        """
        **Inputs**
        
        new_par_dict : dictionary
            Dictionary keys are the parameter names
            Dictionary values are test values for each parameter

        **Returns**

        True if all of the test parameter values are within the limits
        set by self.lims
        """
        for k in new_par_dict.keys():
            if self._check_lims(new_par_dict[k], self.lims[k]):
                pass
            else:
                return False
        return True

    def __getitem__(self, key):
        """
        **Input**

        key : string

        **Returns**

        Parameter value matching `key` as an Astropy.units.Quantity
        """
        assert key in self.keys, print("{} not a valid parameter key".format(key))
        return u.Quantity(self.vals[key], self.units[key])

    def update(self, new_dict):
        """
        **Input**

        new_vals : dict
            Key-value pairs for parameters that you wish to update

        Updates the model parameter values (Model.par_vals)

        If the new value is beyond the model limits, that value is set
        to the nearest limit.
        """
        # Finds the nearest limit
        def find_nearest_lim(v, l):
            if (v < l[0]): return l[0]
            if (v > l[1]): return l[1]
        # Set the new parameters
        for k in new_dict.keys():
            v = new_dict[k]
            if self._check_lims(v, self.lims[k]):
                self.vals[k] = v
            else:
                self.vals[k] = find_nearest_lim(v, self.lims[k])
        return

    # Place-holder function
    def calculate(self, ener_lo, ener_hi):
        """
        Placeholder function. Returns None.
        """
        return

    # Print information about this model
    def info(self):
        """
        Prints model parameters, values, limits, and units.
        """
        print("\n" + "-" * 80)
        print("{}".format(self.name))
        print("{:15}{:10}{:23}{:10}".format('Parameter','Value','Limits','Unit'))
        print("-" * 80)
        for k in self.keys:
            print("{:10}{:10}{:10}{:10}\t{:10}".format(
            k, self.vals[k], self.lims[k][0], self.lims[k][0], self.units[k]
            ))
        return
