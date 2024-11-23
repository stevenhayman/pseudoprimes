"""Module for primality testing and factorization using trial division.

This module provides functions for determining whether a number is prime and
finding its prime factorization using the trial division method.
"""
from math import isqrt
from typing import Dict, Tuple


def is_prime(n: int) -> bool:
    """Test if a number is prime using trial division.

    Args:
        n: The integer to test for primality

    Returns:
        bool: True if n is prime, False otherwise

    Example:
        >>> is_prime(7)
        True
        >>> is_prime(8)
        False
    """
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    m: int = isqrt(n) + 1
    for d in range(3, m, 2):
        if n % d == 0:
            return False

    return True


def factor(n: int) -> Dict[int, int]:
    """Compute the prime factorization of a positive integer.

    Uses trial division to find prime factors. Returns a dictionary where
    keys are prime factors and values are their exponents.

    Args:
        n: The positive integer to factorize

    Returns:
        Dict[int, int]: Prime factors mapped to their exponents

    Raises:
        ValueError: If n is less than 1

    Example:
        >>> factor(12)
        {2: 2, 3: 1}  # represents 2^2 * 3^1 = 12
    """
    if n < 1:
        raise ValueError(f"{n} is not a positive integer")

    factors = {}
    d = 2
    while n % d == 0:
        (n, d, e) = __get_factor(n, d)
        factors[d] = e

    d = 3
    while d ** 2 <= n and n > 1:
        if n % d == 0:
            (n, d, e) = __get_factor(n, d)
            factors[d] = e
        d += 2

    if n != 1:
        factors[n] = 1

    return factors


def __get_factor(n: int, d: int) -> Tuple[int, int, int]:
    """Helper function to extract a prime factor and its exponent.

    Args:
        n: The number being factored
        d: The potential prime factor to test

    Returns
        Tuple[int, int, int]: A tuple containing the remaining number, the prime factor, and the exponent
    """
    e: int = 0
    while n % d == 0:
        n //= d
        e = e + 1
    return (n, d, e)
