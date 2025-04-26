"""This implements simulations algorithms nMGA and LGA, as a Python translation of the implmentation of the work done by Masuda and Rocha in their paper outlinging LGA.

Here we implement using some gillespymax internal tools in order to make comparisons relatively fair.
We note that their algorithm focus on the relative expense of geenrating new random variates, especially since they
use the Mersenne Twister + a floating point division to acquire their Uniform(0,1) random variate.
Relatively speaking, random variate generation is not as expensive in Python, since the (non-random number generation) operations are slower by orders of magnitude.
"""

import numpy as np
import networkx as nx
from typing import Mapping, Hashable, Iterable
from abc import abstractmethod

from .sim import Simulator


class NGMA(Simulator):
    def __init__(
        self,
        graph: nx.Graph,
        initial_state: Mapping[Hashable, str],
        initial_time=0,
        parameters: Mapping | None = None,
        return_statuses: Iterable | None = None,
    ):
        super().__init__(
            graph=graph,
            initial_time=initial_time,
            initial_state=initial_state,
            parameters=parameters,
            return_statuses=return_statuses,
        )

        self.dynamic_rates = dict()


    def update_dynamic_rates(self):
        # extract time, compute the dynamic rates 

class LaplaceGillepsie(Simulator):
    def __init__(
        self,
        graph: nx.Graph,
        initial_state: Mapping[Hashable, str],
        initial_time = 0,
        parameters: Mapping | None = None,
        return_statuses: Iterable | None = None,
    ):
        super().__init__(
            graph=graph,
            initial_time=initial_time,
            initial_state=initial_state,
            parameters=parameters,
            return_statuses=return_statuses,
        )

    @abstractmethod
    def draw_process_rate(self, node, event_class):
        """Draw a new process rate for a given node"""
