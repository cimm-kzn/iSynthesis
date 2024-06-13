# -*- coding: utf-8 -*-
#
#  Copyright 2019-2023 Adelia Fatykhova <adelik21979@gmail.com>
#  This file is part of iSynthesis.
#
#  iSynthesis is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.


from CGRtools.containers import ReactionContainer, MoleculeContainer
from collections import Counter
from iSynthesis.config import db
from hashlib import md5
from pandas import concat
from pony.orm import db_session
from pickle import load, dump
from re import search
from traceback import format_exc

from CGRtools._functions import tuple_hash
from collections import Counter, defaultdict
from StructureFingerprint import LinearFingerprint
from itertools import zip_longest
from .load_data.load_index import index
from math import log2
import faiss
import numpy as np


def query2molecule(queries):
    molecules = []
    for q in queries:
        m = MoleculeContainer()
        for a in q.atoms():
            m.add_atom(a[1].atomic_symbol, a[0])
        for b in q.bonds():
            # m.add_bond(*b)
            m.add_bond(b[0], b[1], b[2].order[0])  #для версии 4.1
        molecules.append(m)
    return molecules

class MyFingerprint(LinearFingerprint):
    def __init__(self,  min_radius: int = 1, max_radius: int = 6, length: int = 4096,
                 number_active_bits: int = 2, number_bit_pairs: int = 4):
        super().__init__(min_radius, max_radius, length, number_active_bits, number_bit_pairs)

    def difference_fingerprint(self, target, reactant, rule):
        rule_dict_reactants1 = {k: v for k, v in self._fragments(query2molecule(rule.reactants)[0]).items()}
        rule_dict_reactants2 = {k: v for k, v in self._fragments(query2molecule(rule.reactants)[1]).items()}

        # reactant1 + reactant2
        for k, v in rule_dict_reactants2.items():
            if k not in rule_dict_reactants1:
                rule_dict_reactants1[k] = v
                continue
            rule_dict_reactants1[k] = v + rule_dict_reactants1[k]

        rule_dict_products = {k: v for k, v in self._fragments(query2molecule(rule.products)[0]).items()}

        # products - reactants
        for k, v in rule_dict_reactants1.items():
            if k not in rule_dict_products:
                rule_dict_products[k] = -v
                continue
            rule_dict_products[k] = rule_dict_products[k] - v

        target_dict = {k: v for k, v in self._fragments(target).items()}
        reactant_dict = {k: v for k, v in self._fragments(reactant).items()}

        for k, v in rule_dict_products.items():
            if k not in reactant_dict:
                
                reactant_dict[k] = v
                continue
            reactant_dict[k] = v + reactant_dict[k]

        for k, v in reactant_dict.items():
            if k not in target_dict:
                target_dict[k] = -v
                continue
            target_dict[k] = target_dict[k] - v
        
        hashes = list({tuple_hash((*tpl, cnt)) for tpl, count in target_dict.items()
                                for cnt in range(min(count, self.number_bit_pairs))})

        number_active_bits = self.number_active_bits
        mask = self.length - 1
        log = int(log2(self.length))
        active_bits = set()

        for tpl in hashes:
            active_bits.add(tpl & mask)
            if number_active_bits == 2:
                active_bits.add(tpl >> log & mask)
            elif number_active_bits > 2:
                for _ in range(1, number_active_bits):
                    tpl >>= log  # shift
                    active_bits.add(tpl & mask)
        
        return active_bits



def find_by_fingerprint(found_fp, operator='substructure'):
    """
    returns ordered list of tuples of all substructure molecules and tanimoto

    """
    if found_fp:
        fp = np.zeros(4096, dtype=bool)
        fp[list(found_fp)] = 1
        fp = np.array(fp)
        query = []
        query.append(fp)
        query = np.array(query).astype('float32')
        
        faiss.normalize_L2(query)
        index.nprobe = 24
        D, I = index.search(query, 100)
        for i, d in zip(I, D):
            res = [k for k in zip(i, d)]
        return res
    return 


def get_reactions(groups, single=True):
    rules_dict = __single if single else __two
    return (r for r, gl in rules_dict.items() if any(g in gl for g in groups))


def index_structure(structure):
    return [i for i, g in __groups.items() if g <= structure]


__single = load(open('iSynthesis/data/rules/single.pickle', 'rb'))
__two = load(open('iSynthesis/data/rules/double.pickle', 'rb'))
__groups = load(open('iSynthesis/data/rules/groups.pickle', 'rb'))


__all__ = ['get_reactions', 'index_structure', 'MyFingerprint', 'find_by_fingerprint']
