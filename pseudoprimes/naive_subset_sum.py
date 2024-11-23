"""Module providing functions for solving the subset sum problem.

This module implements a naive approach to finding subsets of elements that sum
to a target value, supporting generic types through callable operations.
"""
from __future__ import annotations

from typing import Callable, Iterable, Optional, Set, TypeVar

from pseudoprimes.sets import powerset

T = TypeVar('T')


def find_subsets(
    elements: Set[T],
    target: T,
    addition: Callable[[T, T], T],
    equals: Callable[[T, T], bool],
    identity_element: T
) -> Iterable[Set[T]]:
    """Find all subsets of elements that sum to the target value.
    
    Args:
        elements: Set of elements to search through
        target: The target sum to find
        addition: Function that defines how to add two elements
        equals: Function that defines equality between two elements
        identity_element: The additive identity element (e.g., 0 for integers)
        
    Returns:
        Iterable[Set[T]]: All subsets whose elements sum to the target value
        
    Example:
        >>> list(find_subsets({1, 2, 3}, 3, int.__add__, int.__eq__, 0))
        [{3}, {1, 2}]
    """
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
    """Find all non-empty subsets of elements that sum to the target value.
    
    Args:
        elements: Set of elements to search through
        target: The target sum to find
        addition: Function that defines how to add two elements
        equals: Function that defines equality between two elements
        identity_element: The additive identity element (e.g., 0 for integers)
        
    Returns:
        Iterable[Set[T]]: All non-empty subsets whose elements sum to the target value
        
    Example:
        >>> list(find_non_empty_subsets({1, 2, 3}, 3, int.__add__, int.__eq__, 0))
        [{3}, {1, 2}]
    """
    for subset in find_subsets(elements, target, addition, equals, identity_element):
        if len(subset) > 0:
            yield subset


def find_first_non_empty_subset(
    elements: Set[T],
    target: T,
    addition: Callable[[T, T], T],
    equals: Callable[[T, T], bool],
    identity_element: T
) -> Optional[Set[T]]:
    """Find the first non-empty subset that sums to the target value.
    
    Args:
        elements: Set of elements to search through
        target: The target sum to find
        addition: Function that defines how to add two elements
        equals: Function that defines equality between two elements
        identity_element: The additive identity element (e.g., 0 for integers)
        
    Returns:
        Optional[Set[T]]: The first non-empty subset found that sums to the target value,
            or None if no such subset exists
        
    Example:
        >>> find_first_non_empty_subset({1, 2, 3}, 3, int.__add__, int.__eq__, 0)
        {3}  # or {1, 2}, depending on iteration order
    """
    for subset in find_non_empty_subsets(elements, target, addition, equals, identity_element):
        return subset

    return None
