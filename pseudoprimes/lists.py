"""Functions for splitting lists into sublists.

This module provides utilities for dividing lists into smaller sublists,
particularly useful for parallel processing and data partitioning.
"""
from typing import List, TypeVar

T = TypeVar('T')


def split_equally(values: List[T], sublists_count: int) -> List[List[T]]:
    """Split a list into approximately equal-sized sublists.
    
    Args:
        values: The input list to be split
        sublists_count: The number of sublists to create
        
    Returns:
        List[List[T]]: A list containing sublists_count sublists, where the sizes
        differ by at most 1 to account for when len(values) is not evenly divisible
        by sublists_count
        
    Example:
        >>> split_equally([1, 2, 3, 4, 5], 2)
        [[1, 2, 3], [4, 5]]
    """
    base_size, remainder = divmod(len(values), sublists_count)

    return [
        values[i * base_size + min(i, remainder): (i + 1)
               * base_size + min(i + 1, remainder)]
        for i in range(sublists_count)
    ]
