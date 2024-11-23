"""Functions for constructing higher-order Carmichael numbers.

This module implements algorithms for generating prime numbers that can be used
to construct higher-order Carmichael numbers based on Erdős's construction method."""
from typing import Iterable, Set
from pseudoprimes.factorization import Factorization
from pseudoprimes.integers import int_root
from pseudoprimes.partitioned_divisors_iterator import PartitionedDivisorsIterator
from pseudoprimes.fellows_koblitz import is_prime


def build_primes_set(factorization: Factorization, order: int) -> Set[int]:
    """Generate a set of primes for constructing a higher-order Carmichael number.
    
    Args:
        factorization: The base factorization to build from
        order: The order of the Carmichael number to construct"""
    limit = int_root(int(factorization), order)
    divisors_iterator = map(
        int,
        _generate_divisors_less_than(factorization, limit)
    )
    return _build_primes_set_for_divisors(
        factorization,
        order,
        divisors_iterator
    )


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


def _is_valid_prime(
    p: int,
    divisor: int,
    factorization: Factorization,
    order: int,
    value: int
) -> bool:
    """Checks if a prime p is valid for Carmichael number construction."""
    return (
        not factorization.has_prime_divisor(p)
        and all(value % (p**r - 1) == 0 for r in range(1, order + 1))
        and is_prime(p, divisor)
    )


def _generate_divisors_less_than(
    factorization: Factorization,
    limit: int
) -> Iterable[Factorization]:
    """Generates divisors of a factorization less than a given limit."""
    if len(factorization.factors()) == 0:
        if int(factorization) < limit:
            yield factorization
    else:
        primes = factorization.primes()
        exponents = factorization.exponents()

        prime = primes[0]
        exponent = exponents[0]

        for i in range(1 if prime == 2 else 0, exponent + 1):
            partial_value = prime ** i
            if partial_value >= limit:
                break
            partial_factors = dict(zip(primes[1:], exponents[1:]))
            partial_factorization = Factorization(factors=partial_factors)

            for partial_divisor in _generate_divisors_less_than(
                partial_factorization,
                limit // partial_value
            ):
                new_primes = [prime] + partial_divisor.primes()
                new_exponents = [i] + partial_divisor.exponents()
                divisor_factors = dict(zip(new_primes, new_exponents))
                divisor = Factorization(factors=divisor_factors)
                yield divisor
