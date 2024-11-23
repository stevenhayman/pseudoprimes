"""Functions for working with primorials.

This module provides utilities for calculating primorials - numbers that are
products of consecutive prime numbers starting from 2.
"""
from pseudoprimes.factorization import Factorization
from pseudoprimes.sieve_of_eratosthenes import primes_less_than


def find_largest_primorial_less_than(n: int) -> Factorization:
    """Find the factorization of the largest primorial less than n.
    
    A primorial is the product of consecutive prime numbers starting from 2.
    For example, the first few primorials are:
    2, 2×3, 2×3×5, 2×3×5×7, ...
    
    Args:
        n: The upper bound to search below
        
    Returns:
        Factorization: The prime factorization of the largest primorial less than n
        
    Example:
        >>> find_largest_primorial_less_than(100)
        # Returns factorization of 2×3×5 = 30, since 2×3×5×7 = 210 > 100
    """
    primes_limit = n.bit_length() - 1
    primes = primes_less_than(primes_limit)
    primorial = 1
    primorial_primes = []

    for prime in primes:
        primorial *= prime
        if primorial < n:
            primorial_primes.append(prime)

    factors = dict(zip(primorial_primes, [1] * len(primorial_primes)))
    return Factorization(factors=factors)
