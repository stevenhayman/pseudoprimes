"""Module providing a function for performing primality testing via trial division."""
from math import isqrt
from typing import Dict


def is_prime(n: int) -> bool:
    """Function returning True if an integer is prime, and False otherwise."""
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
    """Function returning the prime factorization of a positive integer"""
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


def __get_factor(n, d):
    e: int = 0
    while n % d == 0:
        n //= d
        e = e + 1
    return (n, d, e)
