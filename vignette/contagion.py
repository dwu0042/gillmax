import random
import math
from collections import OrderedDict, defaultdict
from enum import Enum, auto
from warnings import warn
import networkx as nx
from scipy.special import gammaincc as upper_incomplete_gamma

from gillespymax import GillespieMaxSim, BaseEvent, NoEvent

from typing import Mapping, Iterable, Hashable, Any, SupportsFloat, Tuple
from os import PathLike


class SimpleContagionSim(GillespieMaxSim):

    _all_states = OrderedDict((
        ('S', 'Susceptible'),
        ('E', 'Exposed'),
        ('I', 'Infected'),
        ('R', 'Recovered'),
        ('T', 'Tested'),
        ('D', 'Dead'),
        ('0', 'group-type node'),
    ))

    _parameters = OrderedDict((
        ('beta', 'transmissivity parameter'),
        ('prop_time_at_home', 'proportion of time an indivdiual spends at home vs in the community'),
        ('sigma_demographic', 'susceptibility by demographic'),
        # ('sigma_group', 'susceptibility by context'),
        ('incubation_scale', 'scale parameter for incubation period (gamma-dist)'),
        ('incubation_shape', 'shape parameter for incubation period (gamma-dist)'),
        ('alpha_recover', 'rate of recovery (non-death)'),
        ('alpha_mort', 'rate of removal (death)'),
        ('prob_death', 'probability of death'),
        ('p_test_0', 'probabililty of test-seeking while symptomatic'),
        ('kappa', 'test return rate'),
    ))

    _sim_objects = {
        'outcome': 'threshold for predetermined outcomes (death/recovery)',
        'gamma_hazard': 'function that returns the actual gamma-distributed hazard at a given time',
        'entry_time': 'time at which an individual entered their current state',
    }

    class Event(BaseEvent, Enum):

        spontaneous = auto()
        infect = auto()
        seek_test = auto()
        get_test_result = auto()
        removal = auto()

    def __init__(
        self,
        graph: nx.Graph,
        initial_state: Mapping[Hashable, str],
        initial_time=0,
        parameters: Mapping | None = None,
        return_statuses: Iterable | None = None,
    ):
        super().__init__(graph=graph, initial_state=initial_state, initial_time=initial_time, parameters=parameters, return_statuses=return_statuses)

        scale = self.parameters['incubation_scale']
        shape = self.parameters['incubation_shape']
        const = 1./(math.gamma(shape) * scale**shape)
        def gamma_haz(t):
            return const * t**(shape-1) * math.exp(-t/scale) / upper_incomplete_gamma(shape, t/scale)

        self.sim_objects = dict()
        self.sim_objects['gamma_hazard'] = gamma_haz
        self.sim_objects['entry_time'] = defaultdict(float)

    @staticmethod
    def create_initial_state(graph, n_seeds, default_state='S', seed_state='E', obsv_state='0', max_attempts=None):

        initial_state = {node: default_state for node in graph.nodes}

        if max_attempts is None:
            max_attempts = 2*n_seeds

        exposed = set()
        graph_groups = [node for node in graph if graph.nodes[node]['bipartite'] == 0]
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
            initial_state[seed] = obsv_state

        return initial_state

    def rate_function(self, node):
        """Determines the maximal rate of reaction for the given node."""
        state = self.status[node]

        if state in 'SDRT0':
            return 0

        if node not in self.sim_objects['outcome']:
            self.sim_objects['outcome'][node] = {
                'D': random.random(),
                'T': random.random(),
            }

        match state:
            case 'E':
                return 1./self.parameters['incubation_scale']
            case 'I':
                return (self.parameters['beta'] 
                        + self.parameters['alpha_recover']
                        + self.parameters['alpha_mort']
                        + self.parameters['kappa']
                    )
            case _:
                raise RuntimeError(f'Unknown state: {state} of node {node}')

    def transition_choice(self, node):
        """Determines the event that will occur, given a reaction is going to occur for a particular node"""

        state = self.status[node]

        if state == 'I':
            is_test = self.sim_objects['outcome'][node]['D'] < self.parameters['p_test_0']
            is_die = self.sim_objects['outcome'][node]['D'] < self.parameters['prob_death']

        roll = random.random() * self.rates[node]
        if state == 'E':
            if roll < self.sim_objects['gamma_hazard'](self.t - self.sim_objects['entry_time'][node]):
                return self.Event.spontaneous, 'I'
            else:
                return NoEvent.no_event, 
        if state == 'I':
            if roll < self.parameters['alpha_mort']:
                if is_die:
                    return self.Event.removal, 'D'
                else:
                    return NoEvent.no_event, 
            elif roll < (self.parameters['alpha_mort'] + self.parameters['kappa']):
                if is_test:
                    return self.Event.seek_test, 
                else:
                    return NoEvent.no_event, 
            elif roll < (self.parameters['alpha_mort'] + self.parameters['kappa'] + self.parameters['alpha_recover']):
                if not is_die:
                    return self.Event.removal, 'R'
                else:
                    return NoEvent.no_event, 
            else:
                return self.Event.infect, 

        warn(f"Unparsed state {state} of {node}", category=RuntimeWarning)
        return NoEvent.no_event, 

    def manage_event(self, event_type, event_info):
        
        actualised_events = []
        influence_set = []

        if event_type is self.Event.spontaneous:
            node, new_status = event_info
            old_status = self.status[node]
            self.status[node] = new_status
            actualised_events.append({
                't': self.t,
                'enode': node,
                'efrom': old_status,
                'eto': new_status
            })
            influence_set.append(node)
            self.sim_objects['entry_time'][node] = self.t
        elif event_type is self.Event.infect:
            node, = event_info
            node_status = self.status[node]
            context = 'HH' if random.random() < self.parameters['prop_time_at_home'] else 'CC'
            group = random.choice()
            neighbour = random.choice(list(self.graph.neighbours(group)))
