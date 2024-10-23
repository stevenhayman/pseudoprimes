"""Functions for working with primorials."""
from pseudoprimes.factorization import Factorization
from pseudoprimes.sieve_of_eratosthenes import primes_less_than


def find_largest_primorial_less_than(n: int) -> Factorization:
    """Returns the factorization of the largest primorial less than n."""
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
