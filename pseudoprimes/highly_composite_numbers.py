"""Functions for generating factorizations of candidate highly composite numbers."""
from typing import List, Optional

from pseudoprimes.factorization import Factorization
from pseudoprimes.primorial import find_largest_primorial_less_than


def int_log(b: int, n: int) -> int:
    """Returns the floor of log base b of n."""
    i = 0
    value = 1
    while value <= n:
        value *= b
        i += 1
    return i - 1


def __generate_factorizations(
    primes: List[int],
    bound: int,
    previous_exponent: Optional[int] = None
) -> List[Factorization]:
    """Recursively generates factorizations with primes, subject to a bound."""
    if len(primes) == 0:
        return [Factorization(factors={})]

    factorizations = []
    prime = primes[0]

    max_exponent_for_bound = int_log(prime, bound)
    max_exponent = min(previous_exponent, max_exponent_for_bound) if previous_exponent is not None \
        else max_exponent_for_bound

    for exponent in range(max_exponent, -1, -1):
        for partial_factorization in __generate_factorizations(
            primes[1:],
            bound // (prime ** exponent),
            exponent
        ):
            factors = {prime: exponent, **
                       partial_factorization.factors().copy()}
            factorizations.append(Factorization(factors=factors))

    return factorizations


def generate_candidate_highly_composite_numbers_less_than(n: int) -> List[Factorization]:
    """Generates factorizations of candidate highly composite numbers less than n."""
    candidates = []

    primes = find_largest_primorial_less_than(n).primes()

    for factorization in __generate_factorizations(primes, n):
        value = int(factorization)
        if value < n and factorization.has_prime_divisor(3):
            candidates.append(factorization)

    return sorted(candidates)
