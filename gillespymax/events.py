"""Defines the base classes for Events so that NoEvent is correct

Taken from [cobin](https://gitlab.com/cma-public-projects/cobin)
"""

from enum import Enum


class BaseEvent(Enum):
    """Enumeration class that defines an ordering over the elements"""

    def __lt__(self, other):
        if isinstance(other, BaseEvent):
            return self.value < other.value
        else:
            raise TypeError(
                f"{other} is not a valid event type (not derived from BaseEvent)"
            )


class NoEvent(BaseEvent, Enum):
    """A shadow event (no reaction) enum. Used as a sentinel object."""

    no_event = 0
