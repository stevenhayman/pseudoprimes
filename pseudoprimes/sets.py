"""Functions for working with sets.

This module provides utilities for generating mathematical set operations,
including the powerset (set of all subsets) of a given set.
"""
from itertools import chain, combinations
from typing import Set, TypeVar


T = TypeVar('T')


def powerset(elements: Set[T]) -> Set[Set[T]]:
    """Generate the powerset of a set.
    
    The powerset is the set of all possible subsets of a set, including the
    empty set and the set itself.
    
    Args:
        elements: The input set
        
    Returns:
        Set[Set[T]]: The powerset of the input set
        
    Example:
        >>> powerset({1, 2})
        {set(), {1}, {2}, {1, 2}}
    """
    return chain.from_iterable(combinations(elements, r) for r in range(len(elements) + 1))
