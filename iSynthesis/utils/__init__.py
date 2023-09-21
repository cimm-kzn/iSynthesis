# -*- coding: utf-8 -*-
from CGRtools.containers import ReactionContainer, MoleculeContainer
from CIMtools.preprocessing import Fragmentor
from collections import Counter
from iSynthesis.config import db
from hashlib import md5
from pandas import concat
from pony.orm import db_session
from pickle import load
from re import search
from traceback import format_exc


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


def difference_fingerprint(target, reactant, rule):
    fragmentor = Fragmentor(header=False, fragment_type=3,
                            workpath='/tmp', min_length=2,
                            max_length=6, useformalcharge=True)
    reaction = ReactionContainer(query2molecule([x for x in rule.reactants if len(x) > 1]),
                                 query2molecule([x for x in rule.products if len(x) > 1]))
    try:
        reactants = fragmentor.transform(reaction.reactants).sum(axis=0).to_frame().transpose()
        products = fragmentor.transform(reaction.products).sum(axis=0).to_frame().transpose()
    except:
        print(format_exc())
        return
    target_fp = fragmentor.transform([target])
    reactant_fp = fragmentor.transform([reactant])

    z = concat([products, reactants], axis=0, sort=False).fillna(0)  # join
    df = (z.iloc[0] - z.iloc[1]).to_frame().transpose().dropna()  # P - R
    ta = concat([target_fp, reactant_fp], axis=0, sort=False).fillna(0)  # join T - A
    t_a = (ta.iloc[0] - ta.iloc[1]).to_frame().transpose().dropna()  # T - A
    t_a = t_a.loc[:, (t_a != 0).any(axis=0)]
    ddd = concat([t_a, df], axis=0, sort=False).fillna(0)
    ddd = (ddd.iloc[0] - ddd.iloc[1]).to_frame().transpose().dropna()  # T - A - DF
    d = ddd.loc[:, (ddd > 0).any(axis=0)]

    ttt = dict(Counter([int(a) for _, a in target.atoms()]).items())
    rcc = dict(Counter([int(a) for _, a in reactant.atoms()]).items())

    atoms_t_r = {}
    for key in rcc:
        e = ttt.get(key, 0)
        dif = e - rcc[key]
        if dif > 0:
            atoms_t_r[key] = dif
        if e:
            ttt.pop(key)
    for k, v in ttt.items():
        atoms_t_r[k] = v

    bonds = {1: '-', 2: '=', 3: '+', 4: '#'}
    patterns = []
    r_bonds = [(r.atom(x[0]).atomic_symbol, r.atom(x[1]).atomic_symbol, x[2].order)
               for r in reaction.reactants for x in r.bonds()]
    p_bonds = [(p.atom(x[0]).atomic_symbol, p.atom(x[1]).atomic_symbol, x[2].order)
               for p in reaction.products for x in p.bonds()]
    to_delete = set(p_bonds).difference(set(r_bonds))
    for bond in to_delete:
        x, y, b = bond
        patterns.append(f"{x}{bonds[b]}{y}")

    seen = set()
    for pattern in patterns:
        reverse = pattern[::-1]
        if pattern not in seen or reverse not in seen:
            seen.update([pattern, reverse])
            for c in d.columns:
                charge = search(r'\&(.*)\&', c)
                if charge:
                    charge = charge.group(1)
                    formatted = c.replace('&{}&'.format(charge), '')
                else:
                    formatted = c
                if pattern in formatted or reverse in formatted:
                    d = d.drop(c, 1)

    bits_map = {}
    active_bits = set()
    mask = 2 ** 12 - 1

    for k in atoms_t_r:
        active_bits.add(k & mask)
        active_bits.add((k >> 5) & mask)
    for f in d.columns:
        prev = []
        for i in range(1, 4 + 1):
            bs = md5(f'{i}_{f}'.encode()).digest()
            bits_map[(f, i)] = prev = [int.from_bytes(bs[r: r + 2], 'big') & mask
                                       for r in range(0, 2 * 2, 2)] + prev

    for k, v in d.loc[0].items():
        if v:
            active_bits.update(bits_map[(k, v if v < 4 else 4)])

    return active_bits


def find_by_fingerprint(fp, operator='substructure'):
    """
    returns ordered list of tuples of all substructure molecules and tanimoto

    """
    if fp:
        with db_session:
            found = db.Molecule.find_substructure_fingerprint(list(fp))
            if found:
                return found.molecules()
            else:
                found = db.Molecule.find_similar_fingerprint(list(fp))
                if found:
                    return found.molecules()


def get_reactions(groups, single=True):
    rules_dict = __single if single else __two
    return (r for r, gl in rules_dict.items() if any(g in gl for g in groups))


def index_structure(structure):
    return [i for i, g in __groups.items() if g <= structure]


__single = load(open('iSynthesis/data/rules/single.pickle', 'rb'))
__two = load(open('iSynthesis/data/rules/double.pickle', 'rb'))
__groups = load(open('iSynthesis/data/rules/groups.pickle', 'rb'))


__all__ = ['get_reactions', 'index_structure', 'difference_fingerprint', 'find_by_fingerprint']
