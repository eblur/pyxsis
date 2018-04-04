import astropy.units as u

class Model(object):
    def __init__(self, par_keys, par_vals, par_lims, par_units):
        """
        Model superclass
        ----------------
        par_keys : list of parameter names
        par_vals : list of initial parameter values
        par_lims : list of tuples containing parameter limits
        par_units : list of parameter units (strings)
        """
        # Check that the parameter values fall within the limits
        # If they don't, return None
        for v,l in zip(par_vals, par_lims):
            assert self._check_lims(v, l)
            assert l[1] > l[0]
        self.keys = par_keys
        self.vals = dict(zip(par_keys, par_vals))
        self.lims = dict(zip(par_keys, par_lims))
        self.units = dict(zip(par_keys, par_units))
        return

    def _check_lims(self, val, lim):
        """
        Parameters
        ----------

        val : float
        lim : tuple (floats)

        Returns
        -------
        True if `val` is within the boundaries of `lim`
        """
        result = (val >= lim[0]) & (val <= lim[1])
        return result

    def __getitem__(self, key):
        """
        Parameters
        ----------
        key : string

        Returns
        -------
        Parameter value matching `key` as an Astropy.units.Quantity
        """
        assert key in self.keys, print("{} not a valid parameter key".format(key))
        return u.Quantity(self.vals[key], self.units[key])

    def update(self, new_dict):
        """
        Parameters
        ----------
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
