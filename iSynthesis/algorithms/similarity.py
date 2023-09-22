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
