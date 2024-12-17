"""Defines the main simulation logic

Users should implement a subclass of GIllespieMaxSim
"""

import random
from warnings import warn
from abc import ABC, ABCMeta, abstractmethod

import networkx as nx
from sortedcontainers import SortedList

from . import config_loader
from .history import ContagionRecords
from .events import BaseEvent, NoEvent
from .ratedict import RateDict

from typing import Mapping, Iterable, Hashable, Any, SupportsFloat, Tuple
from os import PathLike


class GillespieMaxSim(ABC):

    _states = dict()
    _parameters = dict()
    _sim_objects = dict()

    def __init__(
        self,
        graph: nx.Graph,
        initial_state: Mapping[Hashable, str],
        initial_time=0,
        parameters: Mapping | None = None,
        return_statuses: Iterable | None = None,
    ):

        self.graph = graph

        self.t = initial_time
        self.parameters = dict() if parameters is None else parameters

        self.status = {node: initial_state[node] for node in self.graph.nodes()}

        self.records = ContagionRecords(return_statuses)
        self.records.set_initial_condition(self.t, self.status)

        self.sim_objects = dict()
        self.event_queue = SortedList()

        self.rates = RateDict()

    @abstractmethod
    def maximum_rate(self, node: Hashable) -> SupportsFloat:
        return 0

    @abstractmethod
    def transition_choice(self, node: Hashable) -> Tuple:
        pass

    @abstractmethod
    def manage_event(
        self, event_type: BaseEvent, event_info: Iterable
    ) -> Tuple[Iterable, Iterable]:
        pass

    @classmethod
    def checked_config_load(cls, config_file: PathLike):
        config = config_loader.load(config_file)

        missing_parameters = []
        for parameter in cls._parameters:
            if parameter not in config["parameters"]:
                missing_parameters.append(parameter)

        if len(missing_parameters) > 0:
            raise ValueError(f"Missing parameters in config file: {missing_parameters}")

        return config
    
    @property
    def default_return_states(self):
        return list(self._states.keys())

    def record(self, events: Iterable[Mapping[str, Any]]):
        for event in events:
            self.records.add(**event)

    def update_influence_set(self, influence_set: Iterable[Hashable]):
        for nd in influence_set:
            weight = self.maximum_rate(nd)
            self.rates.insert(nd, weight=weight, cast=float)

    def compute_initial_rates(self):
        for node in self.graph:
            self.rates.insert(node, weight=self.maximum_rate(node))

    def run(self, until=100):

        while self.rates.is_active() or len(self.event_queue):

            # handle event queue
            candidate_delay = self.rates.next_time()
            # extract queued events until none happen before the candidate time to next event
            while len(self.event_queue) > 0 and self.event_queue[0][0] < (
                self.t + candidate_delay
            ):
                t_cand, event_type, *event_info = self.event_queue.pop(0)
                if t_cand < self.t:
                    warn(
                        f"Event Queue produced event in the past ({t_cand}) < ({self.t}, ignoring..."
                    )
                    continue
                self.t = t_cand
                if self.t > until:
                    break
                events, influence_set = self.manage_event(event_type, event_info)
                self.record(events)
                self.update_influence_set(influence_set)
                candidate_delay = self.rates.next_time()

            self.t += candidate_delay

            if self.t >= until:
                break

            # Determine the type of event occuring
            node = self.rates.choose_random()
            event_type, *aux_info = self.transition_choice(node)
            event_info = [node, *aux_info]

            # Manage reaction outcomes
            if event_type is not NoEvent.no_event:
                events, influence_set = self.manage_event(event_type, event_info)
                self.record(events)
                self.update_influence_set(influence_set)

    def write(self, write_to: str | PathLike):
        self.records.write(write_to)
