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


from CGRtools.reactor import Reactor
from logging import getLogger
from .similarity import tversky, tanimoto
from ..utils import difference_fingerprint, find_by_fingerprint
from pony.orm import db_session
from traceback import format_exc


@db_session
def worker(input_queue, output_queue):
    for func, args in iter(input_queue.get, 'STOP'):
        output_queue.put(func(*args))


@db_session
def calc(done_queue, task_number, target):
    logger = getLogger("Synthesis.react.calculate")
    s, to = 0, 0
    for i in range(task_number):
        try:
            res = done_queue.get()
            if res:
                r, t, i = res
                for product in r.products:
                    yield product, tanimoto(product, target), tversky(product, target), r, t
                if i == '1':
                    s += 1
                else:
                    to += 1
        except Exception as e:
            logger.info(f"{e}")
            print(format_exc())
            continue
    logger.info(f"single done: {s}")
    logger.info(f"multi done: {to}")


def react(reactant, template):
    reactor = Reactor(template, delete_atoms=True)
    reaction = next(reactor([reactant]), None)
    if reaction:
        return reaction, template, '1'


def react2mol(target, reactant, template):
    df = difference_fingerprint(target, reactant, template)
    found = find_by_fingerprint(df)
    if found:
        reactor = Reactor(template, delete_atoms=True)
        for i in found:
            with db_session:
                second_reactant = i.structure
                try:
                    reaction = next(reactor([reactant, second_reactant]), None)
                    if reaction:
                        return reaction, template, '2'
                    else:
                        continue
                except Exception as e:
                    logger = getLogger("Synthesis.react2mol")
                    logger.error(e)
