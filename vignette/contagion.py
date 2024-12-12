import random
import math
from collections import OrderedDict, defaultdict
from enum import Enum, auto
from warnings import warn
import networkx as nx
from scipy.special import gammaincc as upper_incomplete_gamma

from gillespymax import GillespieMaxSim, BaseEvent, NoEvent, config_loader

from typing import Mapping, Iterable, Hashable, Any, SupportsFloat, Tuple
from os import PathLike


class SimpleContagionSim(GillespieMaxSim):

    _all_states = OrderedDict(
        (
            ("S", "Susceptible"),
            ("E", "Exposed"),
            ("I", "Infected"),
            ("R", "Recovered"),
            ("T", "Tested"),
            ("Q", "Isolating"),
            ("D", "Dead"),
            ("0", "group-type node"),
        )
    )

    _parameters = OrderedDict(
        (
            ("beta", "transmissivity parameter"),
            (
                "prop_time_at_home",
                "proportion of time an indivdiual spends at home vs in the community",
            ),
            ("sigma_demographic", "susceptibility by demographic"),
            # ('sigma_group', 'susceptibility by context'),
            ("incubation_scale", "scale parameter for incubation period (gamma-dist)"),
            ("incubation_shape", "shape parameter for incubation period (gamma-dist)"),
            ("alpha_recover", "rate of recovery (non-death)"),
            ("alpha_mort", "rate of removal (death)"),
            ("prob_death", "probability of death"),
            ("p_test_0", "probabililty of test-seeking while symptomatic"),
            ("kappa", "test seeking rate"),
            ("test_return", "test return time parameters"),
        )
    )

    _sim_objects = {
        "outcome": "threshold for predetermined outcomes (death/recovery)",
        "gamma_hazard": "function that returns the actual gamma-distributed hazard at a given time",
        "entry_time": "time at which an individual entered their current state",
    }

    class Event(BaseEvent, Enum):

        spontaneous = auto()
        infect = auto()
        seek_test = auto()

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
            initial_state=initial_state,
            initial_time=initial_time,
            parameters=parameters,
            return_statuses=return_statuses,
        )

        scale = self.parameters["incubation_scale"]
        shape = self.parameters["incubation_shape"]
        const = 1.0 / (math.gamma(shape) * scale**shape)

        def gamma_haz(t):
            return (
                const
                * t ** (shape - 1)
                * math.exp(-t / scale)
                / upper_incomplete_gamma(shape, t / scale)
            )

        self.sim_objects = dict()
        self.sim_objects["outcome"] = dict()
        self.sim_objects["gamma_hazard"] = gamma_haz
        self.sim_objects["entry_time"] = defaultdict(float)

        self.compute_initial_rates()

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

    @staticmethod
    def create_initial_state(
        graph,
        n_seeds,
        default_state="S",
        seed_state="E",
        obsv_state="0",
        max_attempts=None,
    ):

        initial_state = {node: default_state for node in graph.nodes}

        if max_attempts is None:
            max_attempts = 2 * n_seeds

        exposed = set()
        graph_groups = [node for node in graph if graph.nodes[node]["bipartite"] == 1]
        for _attempt in range(max_attempts):
            group = random.choice(graph_groups)

            if len(graph[group]) < 1:
                # empty group
                continue

            indv = random.choice(list(graph[group]))
            exposed.add(indv)

            if len(exposed) >= n_seeds:
                break

        for seed in exposed:
            initial_state[seed] = seed_state

        for group in graph_groups:
            initial_state[group] = obsv_state

        return initial_state

    def maximum_rate(self, node):
        """Determines the maximal rate of reaction for the given node."""
        state = self.status[node]

        if state in "SDRTQ0":
            return 0

        if node not in self.sim_objects["outcome"]:
            self.sim_objects["outcome"][node] = {
                "death": random.random(),
                "test_seeking": random.random(),
            }

        match state:
            case "E":
                return 1.0 / self.parameters["incubation_scale"]
            case "I":
                return (
                    self.parameters["beta"]
                    + self.parameters["alpha_recover"]
                    + self.parameters["alpha_mort"]
                    + self.parameters["kappa"]
                )
            case _:
                raise RuntimeError(f"Unknown state: {state} of node {node}")

    def transition_choice(self, node):
        """Determines the event that will occur, given a reaction is going to occur for a particular node"""

        state = self.status[node]

        if state == "I":
            # here, we could save compute by moving the bool computation to rate_function
            # since the quantites are static
            # but this structure allows for cases where the RHS parameters change over time
            is_test = (
                self.sim_objects["outcome"][node]["test_seeking"]
                < self.parameters["p_test_0"]
            )
            is_die = (
                self.sim_objects["outcome"][node]["death"]
                < self.parameters["prob_death"]
            )

        roll = random.random() * self.rates[node]
        if state == "E":
            # rejection sampling / thinning step for non-exponential hazard
            if roll < self.sim_objects["gamma_hazard"](
                self.t - self.sim_objects["entry_time"][node]
            ):
                return self.Event.spontaneous, "I"
            else:
                return (NoEvent.no_event,)
        if state == "I":
            if roll < self.parameters["alpha_mort"]:
                if is_die:
                    return self.Event.spontaneous, "D"
                else:
                    return (NoEvent.no_event,)
            elif roll < (self.parameters["alpha_mort"] + self.parameters["kappa"]):
                if is_test:
                    return (self.Event.seek_test,)
                else:
                    return (NoEvent.no_event,)
            elif roll < (
                self.parameters["alpha_mort"]
                + self.parameters["kappa"]
                + self.parameters["alpha_recover"]
            ):
                if not is_die:
                    return self.Event.spontaneous, "R"
                else:
                    return (NoEvent.no_event,)
            else:
                return (self.Event.infect,)

        warn(f"Unparsed state {state} of {node}", category=RuntimeWarning)
        return (NoEvent.no_event,)

    def change_state(self, node, to_state, **aux_info):
        from_state = self.status[node]
        self.status[node] = to_state
        self.sim_objects["entry_time"][node] = self.t
        return {
            "t": self.t,
            "enode": node,
            "efrom": from_state,
            "eto": to_state,
            **aux_info,
        }

    def manage_event(self, event_type, event_info):

        actualised_events = []
        influence_set = []

        if event_type is self.Event.spontaneous:
            node, new_status = event_info
            state_change = self.change_state(
                node,
                to_state=new_status,
            )
            actualised_events.append(state_change)
            influence_set.append(node)
        elif event_type is self.Event.infect:
            (node,) = event_info
            node_status = self.status[node]
            if node_status == "I":
                context = (
                    "HH"
                    if random.random() < self.parameters["prop_time_at_home"]
                    else "CC"
                )
                all_groups = self.graph.neighbors(node)
                valid_groups = [
                    group for group in all_groups if group.startswith(context)
                ]
                group = random.choice(valid_groups)
                neighbour = random.choice(list(self.graph.neighbors(group)))
                if self.status[neighbour] == "S":
                    # draw for demographic susceptibility
                    neighbour_demography = self.graph.nodes[node].get("demography", 0)
                    if (
                        random.random()
                        < self.parameters["sigma_demographic"][neighbour_demography]
                    ):
                        state_change = self.change_state(
                            neighbour,
                            to_state="E",
                            anode=node,
                            astatus="I",
                            group=group,
                        )
                        actualised_events.append(state_change)
                        influence_set.append(neighbour)
        elif event_type is self.Event.seek_test:
            (node,) = event_info
            if self.status[node] == "I":
                state_change = self.change_state(node, to_state="T")
                self.event_queue.add(
                    (
                        self.t + random.uniform(*self.parameters["test_return"]),
                        self.Event.spontaneous,
                        node,
                        "Q",
                    )
                )
                actualised_events.append(state_change)
                influence_set.append(node)

        return actualised_events, influence_set
