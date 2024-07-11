import random
import math
from collections import OrderedDict, defaultdict
from enum import Enum, auto
from warnings import warn
import networkx as nx
from scipy.special import gammaincc as upper_incomplete_gamma

from networkcontagion.lib.gillespie_max import BaseEvent, NoEvent
from networkcontagion.lib.config_loader import ParameterTable

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
    ('sigma_demographic', 'susceptibility by demographic'),
    ('sigma_group', 'susceptibility by context'),
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
    'entry_time': 'time at which an individual entered their current state',
}

class Event(BaseEvent, Enum):

    spontaneous = auto()
    infect = auto()
    seek_test = auto()
    get_test_result = auto()
    removal = auto()

def rate_function(G, node, status, parameters, sim_objects):
    """Determines the maximal rate of reaction for the given node."""
    state = status[node]

    if state in 'SDRT0':
        return 0

    if node not in sim_objects['outcome']:
        sim_objects['outcome'][node] = {
            'D': random.random(),
            'T': random.random(),
        }

    match state:
        case 'E':
            return 1./parameters['incubation_scale']
        case 'I':
            return (parameters['beta'] 
                    + parameters['alpha_recover']
                    + parameters['alpha_mort']
                    + parameters['kappa']
                )
        case _:
            raise RuntimeError(f'Unknown state: {state} of node {node}')

def transition_choice(G, node, t, status, node_rate, parameters, sim_objects):
    """Determines the event that will occur, given a reaction is going to occur for a particular node"""

    state = status[node]

    if state == 'I':
        is_test = sim_objects['outcome'][node]['D'] < parameters['p_test_0']
        is_die = sim_objects['outcome'][node]['D'] < parameters['prob_death']

    roll = random.random() * node_rate
    if state == 'E':
        if roll < sim_objects['gamma_hazard'](t - sim_objects['entry_time'][node]):
            return Event.infect, 'I'
        else:
            return NoEvent.no_event, 
    if state == 'I':
        if roll < parameters['alpha_mort']:
            if is_die:
                return Event.removal, 'D'
            else:
                return NoEvent.no_event, 
        elif roll < (parameters['alpha_mort'] + parameters['kappa']):
            if is_test:
                return Event.seek_test, 
            else:
                return NoEvent.no_event, 
        elif roll < (parameters['alpha_mort'] + parameters['kappa'] + parameters['alpha_recover']):
            if not is_die:
                return Event.removal, 'R'
            else:
                return NoEvent.no_event, 
        else:
            return Event.infect, 

    warn(f"Unparsed state {state} of {node}", category=RuntimeWarning)
    return NoEvent.no_event, 

def manage_event(G, t, event_type, event_info, status, event_queue, parameters, sim_objects):
    pass



def initialise_parameters(parameters, sim_objects):

    scale = parameters['incubation_scale']
    shape = parameters['incubation_shape']
    const = 1./(math.gamma(shape) * scale**shape)
    def gamma_haz(t):
        return const * t**(shape-1) * math.exp(-t/scale) / upper_incomplete_gamma(shape, t/scale)

    sim_objects['gamma_hazard'] = gamma_haz

    sim_objects['entry_time'] = sim_objects.get('entry_time', defaultdict(float))