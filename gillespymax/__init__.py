"""Gillespie-Max algorithm and library for simulating arbitrary stochastic systems"""

__version__ = "0.0.0a1"

from .sim import GillespieMaxSim
from .events import BaseEvent, NoEvent
from .config_loader import load
