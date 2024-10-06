"""Module providing a function which sieves for primes."""
from math import isqrt
from typing import Set


def primes_less_than(n: int) -> Set[int]:
    """Function returning all primes less than n"""
    is_prime = [True] * n

    for potential_prime in range(2, isqrt(n) + 1):
        if is_prime[potential_prime]:
            for prime_multiple in range(potential_prime ** 2, n, potential_prime):
                is_prime[prime_multiple] = False

    return {potential_prime for potential_prime in range(2, n) if is_prime[potential_prime]}
