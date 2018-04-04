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
        self.keys = par_keys
        self.vals = dict(zip(par_keys, par_vals))
        self.lims = dict(zip(par_keys, par_vals))
        self.units = dict(zip(par_keys, par_units))
        return

    def __getitem__(self, key):
        assert key in self.keys, print("{} not a valid parameter key".format(key))
        return u.Quantity(self.vals[key], self.units[key])
