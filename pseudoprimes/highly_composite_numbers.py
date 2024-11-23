"""Functions for generating highly composite numbers and their factorizations.

This module provides utilities for generating candidate highly composite numbers,
which are positive integers that have more divisors than any smaller positive integer.
The implementation uses prime factorizations to efficiently generate candidates.
"""
from typing import List, Optional

from pseudoprimes.factorization import Factorization
from pseudoprimes.integers import int_log
from pseudoprimes.primorial import find_largest_primorial_less_than


def generate_candidate_highly_composite_numbers_less_than(n: int) -> List[Factorization]:
    """Generate factorizations of candidate highly composite numbers less than n.
    
    Args:
        n: Upper bound for the generated numbers
        
    Returns:
        List[Factorization]: A sorted list of factorizations representing candidate
            highly composite numbers. Only includes numbers that have 3 as a prime
            divisor and are less than n.
    """
    candidates = []

    primes = find_largest_primorial_less_than(n).primes()

    for factorization in __generate_factorizations(primes, n):
        value = int(factorization)
        if value < n and factorization.has_prime_divisor(3):
            candidates.append(factorization)

    return sorted(candidates)


def __generate_factorizations(
    primes: List[int],
    bound: int,
    previous_exponent: Optional[int] = None
) -> List[Factorization]:
    """Recursively generate factorizations using given primes up to a bound.
    
    Args:
        primes: List of prime numbers to use in factorizations
        bound: Upper bound for the generated numbers
        previous_exponent: Optional maximum exponent from previous recursion level
        
    Returns:
        List[Factorization]: List of all valid factorizations that can be constructed
            using the given primes and satisfying the bound constraint
            
    Note:
        This is a helper function that implements the recursive algorithm for
        generating factorizations. Each recursion level handles one prime number
        and tries all valid exponents for that prime.
    """
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
