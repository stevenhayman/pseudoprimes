"""Functions for primality testing using Fellows-Koblitz method.

This module implements primality testing algorithms based on the factorization of n - 1,
following the Fellows-Koblitz approach which examines the multiplicative order of elements
modulo n.
"""
import math
from pseudoprimes.factorization import Factorization, lcm
from pseudoprimes.integers import sub_mod


def is_prime(n: int, factorization: Factorization) -> bool:
    """Test if n is prime using the Fellows-Koblitz primality test.
    
    Args:
        n: The number to test for primality
        factorization: The factorization of n - 1
        
    Returns:
        bool: True if n is prime, False otherwise
        
    Note:
        This implementation uses the multiplicative order of elements modulo n
        and the factorization of n - 1 to determine primality.
    """
    if n in (2, 3):
        return True
    if n < 2 or n % 2 == 0:
        return False

    h = None

    # TODO(steven.paul.hayman@gmail.com): Precision?
    for a in range(2, math.floor(math.log(n) ** 2) + 1):
        if pow(a, n - 1, n) != 1:
            return False
        order_factorization = order(a, n, factorization)
        order_value = int(order_factorization)
        for q in order_factorization.primes():
            if math.gcd(sub_mod(pow(a, order_value // q, n), 1, n), n) > 1:
                return False
        h = order_factorization if h is None else lcm(h, order_factorization)
    return int(h) >= math.isqrt(n) + 1


def order(a: int, n: int, factorization: Factorization) -> Factorization:
    """Calculate the multiplicative order of a modulo n.
    
    Args:
        a: The element whose order we want to compute
        n: The modulus
        factorization: The factorization of n - 1
        
    Returns:
        Factorization: The prime factorization of the multiplicative order of a modulo n
    """
    primes = factorization.primes()
    exponents = factorization.exponents()
    factors = {prime: __get_order_exponent(
        a, n, prime, exponents[i]) for i, prime in enumerate(primes)}
    return Factorization(factors=factors)


def __get_order_exponent(a: int, n: int, prime: int, exponent: int) -> int:
    """Calculate the exponent of a prime in the order factorization.
    
    Args:
        a: The element whose order we're computing
        n: The modulus
        prime: The prime factor being considered
        exponent: The maximum possible exponent for this prime
        
    Returns:
        int: The exact exponent of this prime in the order factorization
    """
    t = pow(a, (n - 1) // (prime ** exponent), n)
    g = 0
    while t != 1:
        t = pow(t, prime, n)
        g += 1
    return g
