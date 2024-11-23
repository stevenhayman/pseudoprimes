"""Module for integer factorization and divisor enumeration.

Provides two main classes:
- Factorization: Represents integers through their prime factorizations
- DivisorsIterator: Iterates through all divisors of a factored integer
"""
from __future__ import annotations

from typing import Any, Dict, List

from pseudoprimes.trial_division import is_prime, factor


class Factorization:
    """Represents an integer through its prime factorization.
    
    Stores factorizations as a dictionary of prime factors and their exponents.
    Example: 12 = 2^2 * 3^1 is stored as {2: 2, 3: 1}
    
    Args:
        value: An integer to be factored, or
        factors: A dict mapping prime factors to exponents
    """

    def __init__(self, **kwargs) -> None:
        if 'value' in kwargs:
            self.__factors = factor(kwargs['value'])
        if 'factors' in kwargs:
            factors = kwargs['factors']
            self.__validate(factors)
            self.__factors = self.__to_standard_form(factors)

    def factors(self) -> Dict[int, int]:
        """Get the prime factorization dictionary.
        
        Returns:
            Dict mapping prime factors to exponents, sorted by prime value
        """
        return self.__factors

    def primes(self) -> List[int]:
        """Get the list of prime factors.
        
        Returns:
            List of prime factors in ascending order
        """
        return list(self.__factors.keys())

    def exponents(self) -> List[int]:
        """Get the list of prime factor exponents.
        
        Returns:
            List of exponents corresponding to each prime factor
        """
        return list(self.__factors.values())

    def has_prime_divisor(self, n: int) -> bool:
        """Check if a number is a prime factor.
        
        Args:
            n: The number to check
            
        Returns:
            True if n is a prime factor in this factorization
        """
        return n in self.__factors

    def __eq__(self, other: Any) -> bool:
        if isinstance(other, Factorization):
            return self.factors() == other.factors()

        return False

    def __int__(self) -> int:
        value: int = 1
        for prime, exponent in self.__factors.items():
            value *= prime ** exponent
        return value

    def __lt__(self, other: Factorization) -> bool:
        return int(self) < int(other)

    def __mul__(self, other: Factorization) -> Factorization:
        """Overloads the multiplication operator for two factorizations."""
        if not isinstance(other, Factorization):
            raise ValueError("Can only multiply with another Factorization")

        factors = self.__factors.copy()

        for prime, exponent in other.factors().items():
            if prime in factors:
                factors[prime] += exponent
            else:
                factors[prime] = exponent

        return Factorization(factors=factors)

    def divisors(self) -> List[Factorization]:
        """Get all positive divisors.
        
        Returns:
            List of all divisors as Factorization objects
        """
        return [divisor for divisor in DivisorsIterator(self)]

    def divisors_count(self) -> int:
        """Calculate the total number of positive divisors.
        
        Returns:
            The count of all positive divisors
        """
        product = 1
        for _, exponent in self.__factors.items():
            product *= (exponent + 1)
        return product

    def __validate(self, factors: Dict[int, int]) -> None:
        for prime, exponent in factors.items():
            if not is_prime(prime):
                raise ValueError(f"{prime} is composite")
            if exponent < 0:
                raise ValueError(f"{exponent} is negative")

    def __to_standard_form(self, factors: Dict[int, int]) -> Dict[int, int]:
        sorted_factors = sorted(factors.items())
        return {prime: exponent for prime, exponent in sorted_factors if exponent > 0}

def lcm(a: Factorization, b: Factorization) -> Factorization:
    """Compute the least common multiple of two factored numbers.
    
    Args:
        a: First factored number
        b: Second factored number
        
    Returns:
        The LCM represented as a Factorization
    """
    factors = {}

    a_factors = a.factors()
    b_factors = b.factors()

    primes = set(a_factors.keys()).union(b_factors.keys())

    for prime in primes:
        exponent = max(a_factors.get(prime, 0), b_factors.get(prime, 0))
        factors[prime] = exponent

    return Factorization(factors=factors)

class DivisorsIterator:
    """Enables iteration through all positive divisors of a factored number.

    Generates divisors by systematically incrementing exponents of prime factors
    up to their maximum values in the original factorization.

    Args:
        factorization (Factorization): The factored number whose divisors to iterate through
    """

    def __init__(self, factorization: Factorization) -> None:
        self.__primes = factorization.primes()
        self.__exponents = factorization.exponents()
        self.__divisor_exponents = None

    def __iter__(self) -> DivisorsIterator:
        return self

    def __next__(self) -> Factorization:
        if self.__divisor_exponents is None:
            self.__divisor_exponents = [0] * len(self.__exponents)
            return Factorization(factors=self.__build_factors())
        for i, _ in enumerate(self.__exponents):
            if self.__divisor_exponents[i] < self.__exponents[i]:
                self.__divisor_exponents[i] += 1
                return Factorization(factors=self.__build_factors())

            self.__divisor_exponents[i] = 0
        raise StopIteration

    def __build_factors(self) -> Dict[int, int]:
        return {self.__primes[i]: self.__divisor_exponents[i] for i in range(len(self.__exponents))}
