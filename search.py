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


from CGRtools.files import SDFRead, SDFWrite
from iSynthesis.algorithms.mcts import MonteCarlo
from iSynthesis.algorithms.preparer import preloaded_tversky, get_reagents
from iSynthesis.config import DEEP, THREADS
from pickle import dump
import argparse
import logging
import os
import subprocess
from datetime import date


current_date = date.today()

parser = argparse.ArgumentParser()
parser.add_argument("target")
parser.add_argument("-s", "--steps", default=10000)
parser.add_argument("-r", "--rnum", default=1000)
parser.add_argument("-o", "--output_dir", default='calcs')
parser.add_argument("--reagents", default="tversky")
parser.add_argument("--cpu", default=12)


args = parser.parse_args()
target = args.target
target_name = target.split('/')[-1].split('.')[0]
print(target_name)
reagents_selection = args.reagents

steps = int(args.steps)
r_num = int(args.rnum)
cpu = int(args.cpu)

logger = logging.getLogger("Synthesis")
logger.setLevel(logging.INFO)
name = f'{target_name}_{current_date}'
output_path = f'{args.output_dir}/{name}'


if not os.path.exists(output_path):
    subprocess.call(f'mkdir -p {output_path}'.split())
fh = logging.FileHandler(f"{output_path}/{name}.log", mode='w')
fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(fh)

target = next(SDFRead(target, ignore_stereo=True))
target.canonicalize()

logger.info(f"Program started. Target: {target_name} {str(target)} steps {steps} reagents {r_num} CPU {cpu}")

g = MonteCarlo(target, target_name, steps, max_depth=5, cpu=cpu, output_dir=output_path)

reagents = preloaded_tversky(target, r_num) if reagents_selection == 'tversky' else get_reagents(target, r_num)

if len(reagents) == 0:
    raise Exception('Reagents not found')

g.start('root', reagents)
g.search()

with open(f'{g._backup_file_name_}', 'wb') as f:
    dump(g, f)
with open(g._file_structures_, 'wb') as s:  # dump structure cache
    dump(g.node_attr_dict_factory._structure_, s)

logger.info("Done! Dumped")
