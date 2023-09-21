# -*- coding: utf-8 -*-
from iSynthesis.config import db, FPS
from iSynthesis.algorithms.similarity import tversky
from iSynthesis.utils import *

from pickle import load
from pony.orm import db_session
from StructureFingerprint import LinearFingerprint


@db_session
def get_reagents(target, number):
    r = db.Molecule.find_similar(target)
    f = r.molecules()
    t = r.tanimotos()
    return [(m.structure, t) for m, t in zip(f, t)][:number]


def preloaded_tversky(target, limit=300, source=None):
    f = load(open(FPS, 'rb'))
    tfp = set(LinearFingerprint(max_radius=6, length=4096, number_bit_pairs=4).transform_bitset([target])[0])
    r = {}
    for mis, mfp in f.items():
        t = tversky(mfp, tfp)
        if t > 0.6:
            r[mis] = t
    del f
    res = {}

    for mis, t in sorted(r.items(), key=lambda x: x[1], reverse=True)[:limit + 1]:
        with db_session:
            res[db.Molecule[mis].structure] = t

    del r
    return {k: v for k, v in sorted(res.items(), key=lambda x: x[1], reverse=True)}


def get_func_groups(reactant):
    # return index_structure(reactant)
    with db_session:
        structure = db.Molecule.find_structure(reactant)
        return [x.id for x in structure.classes] if structure else index_structure(reactant)


__all__ = ['preloaded_tversky', 'get_func_groups']
