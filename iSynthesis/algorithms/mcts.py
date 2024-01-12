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


from datetime import date
from iSynthesis.config import ROOT, THREADS
from iSynthesis.utils import *
from logging import getLogger
from multiprocessing import Process, Queue
from networkx import shortest_path
from operator import itemgetter
from pickle import load, dump, _Pickler
from .similarity import tversky
from .preparer import get_func_groups
from pony.orm import db_session
from .react import calc, react, react2mol, worker

from .tree import Tree


class MonteCarlo(Tree, _Pickler):

    def __init__(self, target, target_name, number=5000, max_depth=5, cpu=THREADS, output_dir='.'):
        """
        :param target: target molecule
        :param number: number of iterations on the Tree
        """
        self.target = target
        self.target_name = target_name
        self.number = number
        self._achieved = {}
        self.len = 100
        self.bb = False
        self.iteration = 1
        self.date = date.today()
        self._file_with_target_name_ = f"{output_dir}/{self.target_name}_with_target{self.date}"
        self._backup_file_name_ = f"{output_dir}/{self.target_name}_backup{self.date}"
        self._file_structures_ = f"{output_dir}/{self.target_name}_structures{self.date}"
        self._file_paths_ = f"{output_dir}/{self.target_name}_paths_{self.date}"
        super().__init__()
        self.max_depth = max_depth
        self.__cpu = cpu

    def __setstate__(self, d):
        self.__dict__.update(d)
        structures = load(open(self._file_structures_, 'rb'))
        self.node_attr_dict_factory._structure_ = structures

    def __getstate__(self):
        return self.__dict__

    @property
    def achieved(self):
        return dict(sorted(self._achieved.items(), key=lambda x: max(x['tanimoto'] for x in x[1]),
                           reverse=True)) if self._achieved else {}

    @achieved.setter
    def achieved(self, node):
        p, tan, r, template = node
        if len(self._achieved) < self.len:
            self._achieved.setdefault(p, []).append({'tanimoto': tan, 'reaction': r, 'template': template})
        else:
            best = list(self.achieved)[0]  # p
            lst = self._achieved[best]  # [{t, r}, ...]
            if any(x['tanimoto'] <= tan for x in lst):
                g = self._achieved.get(p)
                if not g:
                    self._achieved.pop(list(self.achieved)[-1])
                self._achieved.setdefault(p, []).append({'tanimoto': tan, 'reaction': r, 'template': template})

    def search(self):
        if self.bb:
            return
        logger = getLogger('Synthesis.MCTS.search')
        task_queue, done_queue = Queue(), Queue()
        procs = [Process(target=worker, args=(task_queue, done_queue))
                 for _ in range(self.__cpu)]
        [p.start() for p in procs]
        for i in range(self.number):
            i += 1
            self.iteration += 1
            logger.info('*** step {} started ***'.format(i))
            print('step {} started'.format(i))
            best = self.select(ROOT)
            if best is None:
                break
            structure = self.data[best]
            b = self._node[best]
            logger.info('best selected: {},  {} , tversky: {}, score: {}'.format(best, structure, b['value'], b['score']))
            if structure != ROOT:
                groups_in_query = get_func_groups(structure)
                one_component = list(get_reactions(groups_in_query))
                two_component = list(get_reactions(groups_in_query, single=False))
                logger.info(f"found 1-c. rules: {len(one_component)}")
                logger.info(f"found 2-c. rules: {len(two_component)}")
                [task_queue.put((react, (structure, template))) for template in one_component]
                [task_queue.put((react2mol, (self.target, structure, template))) for template in two_component]
                calculated = calc(done_queue, len(one_component) + len(two_component), self.target)
                new_nodes = sorted(set(calculated), key=itemgetter(1), reverse=True)
                done = self.expansion(best, new_nodes)
                self.update_achieved(new_nodes)
                self.backup(best)
                if done:
                    with open(self._file_with_target_name_, 'wb') as f:
                        dump(self, f)
                    with open(self._file_structures_, 'wb') as s:  # dump structure cache
                        dump(self.node_attr_dict_factory._structure_, s)
                    logger.info('dumped into {}'.format(self._file_with_target_name_))
            for x in list(self.achieved.items()):
                logger.info('the best yet: {} with {}\n\n'.format(x[0], x[1][0]['tanimoto']))
                break
            if i % 100 == 0:
                with open(self._backup_file_name_, 'wb') as f:
                    dump(self, f)
                with open(self._file_structures_, 'wb') as s:  # dump structure cache
                    dump(self.node_attr_dict_factory._structure_, s)
                logger.info('SIMILAR PATHS')
                self.show_similar_paths(100)
                logger.info('------------------------------------------')
                logger.info('backup dumped into {}'.format(self._backup_file_name_))
        [p.terminate() for p in procs]
        return self

    def start(self, root, nodes):
        logger = getLogger('Synthesis.MCTS.start')
        self.add_node(root, data=root, score=0, value=0.01, parent=None, visits=0)
        for node in nodes:
            num = len(self.nodes) + 1
            predicted = tversky(node, self.target)
            if node == self.target:
                logger.info('target already exists in the database of building blocks')
                print('target already exists in the database of building blocks', 'magenta')
                self.bb = True
                return
            mean = self.set_mean_value(predicted, 1)
            self.add_node(num, data=node, value=predicted, score=self.set_score(mean, 1),
                          parent=0, visits=0, depth=0, mean_value=mean)
            self.add_edge(root, num)
        return self

    def update_achieved(self, nodes):
        for i in nodes:
            p, tan, _, r, template = i
            exist_templates = [x['template'] for x in self.achieved.get(p, [])]
            exist_reactions = [x['reaction'] for x in self.achieved.get(p, [])]
            if template not in exist_templates or r not in exist_reactions:
                self.achieved = (p, tan, r, template)

    def similar_path(self, number):
        crop = {k: self.achieved[k] for k in list(self.achieved)[:number]}.items()
        return [(shortest_path(self, ROOT, node), tanimoto, reaction)
                for product, list_info in crop for inf in list_info
                for node in self.get_nodes(product)
                for tanimoto, reaction in zip([inf['tanimoto']], [inf['reaction']])]

    def show_similar_paths(self, number):
        """
        :param number: number of similar molecules to show path
        :return: dictionary of similar to target molecules with Tv index and all reactions with templates
        {(mol, Tv): (*reactions, *t}emplates)}

        Example:
        1.0 - Tanimoto
        COc1ccc(CCNC(C)=O)cc1OC>>c1c(OC)c(ccc1CCN)OC [rule]
        c1c(OC)c(ccc1CCN)OC.c1ccc(OCC2OC2)cc1C>>c1c(OCC(O)CNCCc2cc(OC)c(cc2)OC)cccc1C [rule]
        $$$ - separator for different pathways of one product
        COc1ccc(cc1OC)CCNC=O>>c1c(OC)c(ccc1CCN)OC [rule]
        c1c(OC)c(ccc1CCN)OC.c1ccc(OCC2OC2)cc1C>>c1c(OCC(O)CNCCc2cc(OC)c(cc2)OC)cccc1C [rule]
        $$$
        END - separator for different products

        0.9928057553956835
        O(c1c(ccc(CC)c1)OCC2OC2)C>>O(c1cc(C)ccc1OCC2OC2)C [rule]
        O(c1cc(C)ccc1OCC2OC2)C.c1c(C(=O)CN)cc(cc1)OC>>C(CNCC(O)COc1c(OC)cc(C)cc1)(c2cccc(c2)OC)=O [rule]
        C(CNCC(O)COc1c(OC)cc(C)cc1)(c2cccc(c2)OC)=O>>c1(OC)c(ccc(C)c1)OCC(O)CNCCc2cccc(c2)OC [rule]
        $$$
        END
        """
        path = self.similar_path(number)
        output = {}
        for t in path:
            tmp_result = []
            ls = [x for x in t[0] if not isinstance(x, list)]
            for i in zip(ls[1:], ls[2:]):
                data = (*self.get_edge_data(*i)['reaction'], *self.get_edge_data(*i)['template'])
                tmp_result.append(data)
            if tmp_result:
                output.setdefault((max(tmp_result[-1][0].products), t[1]), []).append(tmp_result)

        with open(self._file_paths_, 'w') as log:
            for mol, v in output.items():
                print(mol[1])
                seen = set()
                for p in v:
                    n = 0
                    for i in p:
                        if n == 0 and bytes(i[0]) in seen:
                            break
                        seen.add(bytes(i[0]))
                        print(*i)
                        log.write(*i)
                        log.write('\n')
                        n += 1
                    else:
                        print('$$$')
                        log.write('$$$')
                print('END\n')
                log.write('END\n')
        return output


__all__ = ['MonteCarlo']
