""" Mapping for a rate dictionary that allows for weighted random selection via cdf sampling 

Taken from [cobin](https://gitlab.com/cma-public-projects/cobin)
"""

from collections import defaultdict
import random


class RateDict(object):
    """
    Adapted from EoN simulation _ListDict_
    Uses a cdf sampling on the weights, then a uniform draw over the nodes
    associated with those weights
    Efficient only if there are a small set of unique weights possible
    """

    def __init__(self):
        self.weights = defaultdict(list)
        self.itemmap = dict()
        self.pdf = defaultdict(int)
        self.total_weight = 0.0

    def __str__(self):
        return f"RateDict[items = {len(self)}, total_weight = {self.total_weight}]"

    def __len__(self):
        return sum(int(v / k) for k, v in self.pdf.items())

    def __contains__(self, item):
        return item in self.itemmap and self.itemmap[item] != 0

    def __getitem__(self, key):
        return self.itemmap[key]

    def is_active(self):
        r"""
        Returns whether or not there are any more non-negative weights left
        """
        return self.total_weight > (min(self.pdf, default=0) / 2)

    def insert(self, item, weight=None, cast=None):
        r"""
        If not present, then inserts the thing (with weight if appropriate)
        if already there, replaces the weight unless weight is 0

        If weight is 0, then it removes the item and doesn't replace.

        cast is an arbitrary function that is applied on weight, if it is not set to None.

        WARNING:
            replaces weight if already present, does not increment weight.
        """
        self.remove(item)
        if weight != 0:
            if cast is not None:
                weight = cast(weight)
            self.weights[weight].append(item)
            self.itemmap[item] = weight
            self.pdf[weight] += weight
            self.total_weight += weight

    def remove(self, item):
        r"""
        Removes a given item, if it exists.
        """
        if item in self:
            w = self.itemmap[item]
            self.weights[w].remove(item)
            self.pdf[w] -= w
            self.total_weight -= w
            del self.itemmap[item]

    def choose_random(self):
        r"""
        Chooses a random node using reverse CDF mapping.
        """
        w = random.choices(list(self.weights.keys()), weights=self.pdf.values(), k=1)
        return random.choice(self.weights[w[0]])

    def random_removal(self):
        r"""
        Randomly choose and remove a node, which is then returned.
        """
        choice = self.choose_random()
        self.remove(choice)
        return choice

    def next_time(self):
        """Draws the time to the next event (global)

        Returns:
            float: time to next event (global)
        """
        if self.total_weight > 0:
            return random.expovariate(self.total_weight)
        else:
            return float("Inf")
