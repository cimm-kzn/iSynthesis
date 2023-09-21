# -*- coding: utf-8 -*-
from StructureFingerprint import LinearFingerprint


def fingerprints_for_index(query, target):
    if isinstance(query, set):
        query_fp, target_fp = query, target
    else:
        query_fp, target_fp = (set(LinearFingerprint(max_radius=6, length=4096, number_bit_pairs=4).transform_bitset([x])[0])
                               for x in [query, target])
    return len(query_fp), len(target_fp), len(query_fp.intersection(target_fp))


def tanimoto(query, target):
    """
    оценка нод.
    возвращает танимото для пары запрос-результат.
    """
    qc, rc, common = fingerprints_for_index(query, target)
    return common / (qc + rc - common)


def tversky(query, target):
    qc, rc, common = fingerprints_for_index(query, target)
    c = (1 - 10 ** (- qc / 20))
    return common / (0.8 * (qc - common) + 0.2 * (rc - common) + common) * c


__all__ = ['tversky', 'tanimoto']
