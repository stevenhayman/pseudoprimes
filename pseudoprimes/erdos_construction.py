"""Module providing functions for constructing Carmichael numbers."""
from typing import Callable, Set
from pseudoprimes.factorization import Factorization
from pseudoprimes.trial_division import is_prime

def construct_carmichael_number(
    factorization: Factorization,
    find_subset: Callable[[Set[int]], Set[int]]
) -> Set[int]:
    """A function which constructs a Carmichael number using the Erdos construction."""
    primes = build_primes_set(factorization)
    return find_subset(primes)

def build_primes_set(factorization: Factorization) -> Set[int]:
    """A function which builds the set of prime numbers used by the Erdos construction."""
    primes: Set[int] = set()
    for divisor in factorization.divisors():
        p = int(divisor) + 1
        if is_prime(p) and not factorization.has_prime_divisor(p):
            primes.add(p)
    return primes