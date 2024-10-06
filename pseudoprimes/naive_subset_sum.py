"""Module providing classes and functions for finding solutions to the subset sum problem."""
from __future__ import annotations

from itertools import chain, combinations
from typing import Callable, Iterable, Optional, Set, TypeVar

T = TypeVar('T')

def powerset(elements: Set[T]) -> Set[Set[T]]:
    """Function which generates the powerset of a set."""
    return chain.from_iterable(combinations(elements, r) for r in range(len(elements) + 1))

def find_subsets(
    elements: Set[T],
    target: T,
    addition: Callable[[T, T], T],
    equals: Callable[[T, T], bool],
    identity_element: T
) -> Iterable[Set[T]]:
    """Function which returns a subset of elements that sum up to the target value."""
    for subset in powerset(elements):
        subset_sum = identity_element
        for element in subset:
            subset_sum = addition(subset_sum, element)

        if equals(subset_sum, target):
            yield set(subset)

def find_non_empty_subsets(
    elements: Set[T],
    target: T,
    addition: Callable[[T, T], T],
    equals: Callable[[T, T], bool],
    identity_element: T
) -> Iterable[Set[T]]:
    """Function which generates the non-empty subsets of elements that sum up to the target value."""
    for subset in find_subsets(elements, target, addition, equals, identity_element):
        if subset:  # Only yield non-empty subsets
            yield subset

def find_first_non_empty_subset(
    elements: Set[T],
    target: T,
    addition: Callable[[T, T], T],
    equals: Callable[[T, T], bool],
    identity_element: T
) -> Optional[Set[T]]:
    """Returns the first non-empty subset of elements that sums up to the target value."""
    for subset in find_non_empty_subsets(elements, target, addition, equals, identity_element):
        return subset

    return None
