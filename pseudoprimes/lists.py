"""Functions for splitting lists into sublists."""
from typing import List, TypeVar

T = TypeVar('T')


def split_equally(values: List[T], sublists_count: int) -> List[List[T]]:
    """Splits a list into sublists of roughly equal size."""
    base_size, remainder = divmod(len(values), sublists_count)

    return [
        values[i * base_size + min(i, remainder): (i + 1)
               * base_size + min(i + 1, remainder)]
        for i in range(sublists_count)
    ]
