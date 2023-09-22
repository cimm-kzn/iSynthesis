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


from collections import MutableMapping
from CGRtools.containers import MoleculeContainer


class StructureCache(MutableMapping):
    _structure_ = {}

    def __init__(self):
        self._data = {}
        self._str_keys = set()

    def __setitem__(self, key, value):
        if isinstance(value, MoleculeContainer):
            self._str_keys.add(key)
            self._data[key] = str(value)
            if str(value) not in self._structure_:
                self._structure_[str(value)] = value
        else:
            self._data[key] = value

    def __delitem__(self, key):
        del self._data[key]
        if key in self._str_keys:
            self._str_keys.discard(key)

    def __iter__(self):
        return iter(self._data)

    def __getitem__(self, key):
        if key in self._str_keys:
            return self._structure_[self._data[key]]
        return self._data[key]

    def __len__(self):
        return len(self._data)
