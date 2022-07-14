# Copyright (c) 2020-2021 Alec White
# Licensed under the MIT License (see LICENSE for details)
from copy import copy


class Idx(object):
    """
    Index class

    Attributes:
        index (int): Integer labelling an index as unique
        space (str): Name of the index space
        fermion (bool): Index of fermiond space?
    """
    def __init__(self, index, space, fermion=True):
        self.index = index
        self.space = space
        self.fermion = fermion

    def __repr__(self):
        return str(self.index) + "(" + self.space + ")"

    def __hash__(self):
        """Used by e.g. set.intersection"""
        return hash(f"{self.index}({self.space})")

    def __eq__(self, other):
        return self.index == other.index and self.space == other.space

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if self.space < other.space:
            return True
        elif self.space == other.space:
            return self.index < other.index
        else:
            return False

    def __le__(self, other):
        return self < other or self == other

    def __gt__(self, other):
        return not self <= other

    def __ge__(self, other):
        return not self < other


def idx_copy(idx):
    """
    Copy an index by making a deep copy of the index (int) and
    fermion (bool) variables and copying the reference to the space.
    """
    if isinstance(idx, NamedIndex):
        return NamedIndex(copy(idx.index), idx.space, copy(idx.name), bool(idx.fermion))
    return Idx(copy(idx.index), idx.space, bool(idx.fermion))


def is_occupied(idx, occ=None):
    """
    occ: list of possible occupied index spaces
    """
    if occ is None:
        return 'o' in idx.space
    else:
        return idx.space in occ


class NamedIndex(Idx):
    """
    Named index class

    Attributes:
        index (int): Integer labelling an index as unique
        space (str): Name of the index space
        fermion (bool): Index of fermiond space?
        name (str): Name of the symbol
    """
    def __init__(self, index, space, name, fermion=True):
        self.index = index
        self.space = space
        self.name = name
        self.fermion = fermion