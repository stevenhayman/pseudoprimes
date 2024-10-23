"""Functions for constructing higher-order Carmichael numbers."""
from typing import Iterable, Set
from pseudoprimes.factorization import Factorization
from pseudoprimes.partitioned_divisors_iterator import PartitionedDivisorsIterator
from pseudoprimes.fellows_koblitz import is_prime


def _is_valid_prime(
    p: int,
    divisor: int,
    factorization: Factorization,
    order: int,
    value: int
) -> bool:
    """Checks if a prime p is valid for Carmichael number construction."""
    return (
        is_prime(p, divisor)
        and not factorization.has_prime_divisor(p)
        and all(value % (p**r - 1) == 0 for r in range(1, order + 1))
    )


def _build_primes_set_for_divisors(
    factorization: Factorization,
    order: int,
    divisors: Iterable[int],
) -> Set[int]:
    """Builds the set of valid primes from divisors."""
    value = int(factorization)

    primes: Set[int] = set()

    for divisor in divisors:
        p = int(divisor) + 1
        if _is_valid_prime(p, divisor, factorization, order, value):
            primes.add(p)

    return primes


def build_primes_set(factorization: Factorization, order: int) -> Set[int]:
    """Returns a set of primes for constructing a higher-order Carmichael number."""
    return _build_primes_set_for_divisors(factorization, order, factorization.divisors())


def build_primes_set_partition(
    factorization: Factorization,
    order: int,
    partition_index: int,
    partitions_count: int
) -> Set[int]:
    """Returns a partitioned set of primes for constructing a higher-order Carmichael number."""
    divisors_iterator = PartitionedDivisorsIterator(
        factorization,
        partition_index,
        partitions_count
    )

    return _build_primes_set_for_divisors(factorization, order, divisors_iterator)
