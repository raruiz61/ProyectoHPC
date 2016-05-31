# enables import of "epymol"
#
#
# ePyMOL: Closed-source add-ons that are included in incentive versions of PyMOL

def extendapi(name, function=None):
    '''Add as PyMOL command and as function to cmd API'''
    from pymol import cmd as _self
    function = _self.extend(name, function)
    setattr(_self, function.__name__, function)
    return function

try:
    from comm import \
        obscure
except ImportError:
    pass
