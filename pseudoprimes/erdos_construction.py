"""Functions for constructing Carmichael numbers using Erdős's method.

This module implements the Erdős construction method for generating Carmichael numbers
by finding sets of prime numbers that satisfy specific conditions.
"""
from typing import Callable, Set
from pseudoprimes.factorization import Factorization
from pseudoprimes.trial_division import is_prime


def construct_carmichael_number(
    factorization: Factorization,
    find_subset: Callable[[Set[int]], Set[int]]
) -> Set[int]:
    """Construct a Carmichael number using the Erdős construction method.
    
    Args:
        factorization: The base factorization to build from
        find_subset: A function that selects a valid subset of primes
        
    Returns:
        A set of prime numbers whose product forms a Carmichael number
    """
    primes = build_primes_set(factorization)
    return find_subset(primes)


def build_primes_set(factorization: Factorization) -> Set[int]:
    """Build a set of prime numbers that can be used in the Erdős construction.
    
    For each divisor d of the factorization, checks if d+1 is prime and not already
    a divisor of the factorization. These primes are candidates for constructing
    Carmichael numbers.
    
    Args:
        factorization: The base factorization to generate primes from
        
    Returns:
        A set of prime numbers that satisfy the initial conditions for
        Carmichael number construction
    """
    primes: Set[int] = set()
    for divisor in factorization.divisors():
        p = int(divisor) + 1
        if is_prime(p) and not factorization.has_prime_divisor(p):
            primes.add(p)
    return primes
