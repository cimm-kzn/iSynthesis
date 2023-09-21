# -*- coding: utf-8 -*-
from ..config import SEED
from iSynthesis.config import EPSILON, ROOT, VALUE, VISITS, SCORE, DATA, DEPTH, MEAN
from logging import getLogger
from math import log10, log, sqrt
from .cache import StructureCache
from networkx import DiGraph
import random


class Tree(DiGraph):
    """
    base class for MCTS
    """

    def __init__(self):
        self.node_attr_dict_factory = StructureCache
        self.max_depth = 5
        super().__init__()

    def get_score(self, node):
        return self.mean_value[node] + sqrt(2) * sqrt((log(self.iteration) / self.visits[node]))

    def set_score(self, mean, visits):
        return mean + sqrt(2) * sqrt(log(self.iteration) / visits)

    def get_mean_value(self, node):
        return - log10(1 - self.value[node] + 0.01) / self.visits[node]

    @staticmethod
    def set_mean_value(value, visits):
        return - log10(1 - value + 0.01) / visits

    @property
    def score(self):
        return self.nodes(data=SCORE)

    @property
    def value(self):
        return self.nodes(data=VALUE)

    @property
    def mean_value(self):
        return self.nodes(data=MEAN)

    @property
    def visits(self):
        return self.nodes(data=VISITS)

    @property
    def data(self):
        return self.nodes(data=DATA)

    @property
    def depth(self):
        return self.nodes(data=DEPTH)

    def get_nodes(self, node):
        return (x for x, y in self.nodes(data=True) if y['data'] == node)

    def set_node_score(self, node, score):
        self.nodes[node][SCORE] = score

    def set_node_visits(self, node):
        self.nodes[node][VISITS] += 1

    def set_node_depth(self, node, d):
        self.nodes[node][DEPTH] += d

    def set_node_mean(self, node, mean):
        self.nodes[node][MEAN] = mean

    def best_child(self, node):
        candidates = [x for x in self.neighbors(node) if x != self.target and self.depth[x] < self.max_depth]
        if candidates:
            return max(candidates, key=lambda i: self.score[i]) or node
        else:
            raise ValueError('no candidates')
            # logger = getLogger('Synthesis')
            # logger.warning('BEST CHILD IS NONE')
            # return self.random_child(node)

    def select(self, node):
        current = node
        random.seed(SEED)
        while list(self.successors(current)):
            current = self.random_child(current) if random.random() < EPSILON else self.best_child(current)
        return current

    def random_child(self, node):
        random.seed(SEED)
        return random.choice([x for x in self.successors(node) if x != self.target])

    def backup(self, node):
        current = node
        depth = self.depth[current]
        while list(self.predecessors(current)):
            self.set_node_visits(current)
            scores = [self.score[x] for x in self.successors(current)]
            if len(scores) and 0 not in scores:
                score = sum(scores) / len(scores)
            else:
                score = self.get_score(current)
            self.set_node_score(current, score)
            mean_values = [self.mean_value[x] for x in self.successors(current)]
            if len(mean_values) and 0 not in mean_values:
                mean = sum(mean_values) / len(mean_values)
            else:
                mean = self.get_mean_value(current)
            self.set_node_mean(current, mean)
            current = list(self.predecessors(current))[0]
            if current != ROOT:
                depth += 1
                self.set_node_depth(current, depth)
        return self

    def expansion(self, parent, nodes):
        done = False
        logger = getLogger("Synthesis.tree.expansion")
        seen = {}
        for structure, _, tv, reaction, template in nodes:
            sig = bytes(structure)
            if sig not in seen:
                num = len(self.nodes) + 1
                seen[sig] = num
                mean = self.set_mean_value(tv, 1)
                self.add_node(num, data=structure, value=tv, score=self.set_score(mean, 1),
                              visits=0, parent=parent, depth=0, mean_value=mean)
                self.add_edge(parent, num, reaction={reaction}, template={template})
                if structure == self.target:
                    done = True
                    logger.info(f"{'=' * 20}TARGET RECEIVED{'=' * 20}")
                    logger.info(f"{str(reaction)}")
                    print(f"{'=' * 33}TARGET RECEIVED{'=' * 33}")
            else:
                num = seen[sig]
                self.edges[parent, num]['template'].add(template)
                self.edges[parent, num]['reaction'].add(reaction)
        return done

    def best_reagent(self):
        try:
            return max([x for x in self._node.items()
                        if self.has_predecessor(x[0], ROOT) and list(self.successors(x[0]))],
                       key=lambda n: n[1][SCORE])
        except ValueError:
            print('empty tree')

