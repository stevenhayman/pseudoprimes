"""Functions for primality testing based on the factorization of n - 1."""
import math
from pseudoprimes.factorization import Factorization


def __get_order_exponent(a, n, prime, exponent) -> int:
    t = pow(a, (n - 1) // (prime ** exponent), n)
    g = 0
    while t != 1:
        t = pow(t, prime, n)
        g += 1
    return g


def order(a: int, n: int, factorization: Factorization) -> Factorization:
    """Returns the order of a in the multiplicative group mod n."""
    primes = factorization.primes()
    exponents = factorization.exponents()
    factors = {prime: __get_order_exponent(
        a, n, prime, exponents[i]) for i, prime in enumerate(primes)}
    return Factorization(factors=factors)


def sub_mod(a: int, b: int, n: int) -> int:
    """Returns (a - b) mod n."""
    t = (a % n) - (b % n)
    if t < 0:
        t += n
    elif t > n:
        t -= n
    return t


def lcm(a: Factorization, b: Factorization) -> Factorization:
    """Returns the least common multiple of two factorizations."""
    factors = {}

    a_factors = a.factors()
    b_factors = b.factors()

    primes = set(a_factors.keys()).union(b_factors.keys())

    for prime in primes:
        exponent = max(a_factors.get(prime, 0), b_factors.get(prime, 0))
        factors[prime] = exponent

    return Factorization(factors=factors)


def is_prime(n: int, factorization: Factorization) -> bool:
    """Returns whether n is prime using the factorization of n - 1."""
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
