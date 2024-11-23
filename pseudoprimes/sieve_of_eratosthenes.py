"""Module implementing the Sieve of Eratosthenes algorithm.

This module provides an efficient implementation of the classical Sieve of
Eratosthenes algorithm for finding all prime numbers up to a given limit.
"""
from math import isqrt
from typing import Set


def primes_less_than(n: int) -> Set[int]:
    """Find all prime numbers less than n using the Sieve of Eratosthenes.
    
    Args:
        n: Upper bound (exclusive) for finding primes
        
    Returns:
        Set[int]: A set containing all prime numbers less than n
        
    Example:
        >>> primes_less_than(10)
        {2, 3, 5, 7}
    """
    is_prime = [True] * n

    for potential_prime in range(2, isqrt(n) + 1):
        if is_prime[potential_prime]:
            for prime_multiple in range(potential_prime ** 2, n, potential_prime):
                is_prime[prime_multiple] = False

    return {potential_prime for potential_prime in range(2, n) if is_prime[potential_prime]}
