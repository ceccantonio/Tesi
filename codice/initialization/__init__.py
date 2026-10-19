"""Initialization of the lab instruments, one module per instrument (laser.py, oscilloscope.py, motors.py).

    from initialization import init
    laser, scope, motors = init()           # init(laser=False, ...) leaves an instrument out (it is None)
"""
from .laser import Laser
from .motors import Motors
from .oscilloscope import Oscilloscope


def init(laser=True, scope=True, motors=True, n_samples=None):
    """Connect the instruments and return them as (laser, scope, motors); the ones left out are None."""
    return (Laser() if laser else None,
            Oscilloscope(n_samples=n_samples) if scope else None,
            Motors() if motors else None)
