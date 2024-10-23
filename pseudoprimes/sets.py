"""Functions for working with sets."""
from itertools import chain, combinations
from typing import Set, TypeVar


T = TypeVar('T')


def powerset(elements: Set[T]) -> Set[Set[T]]:
    """Generates the powerset of a set."""
    return chain.from_iterable(combinations(elements, r) for r in range(len(elements) + 1))
