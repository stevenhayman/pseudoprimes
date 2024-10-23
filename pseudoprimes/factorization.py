"""Module providing a class for storing the factorization of an integer
and a class for interating through the divisors of a factored integer."""
from __future__ import annotations

from typing import Any, Dict, List

from pseudoprimes.trial_division import is_prime, factor


class Factorization:
    """Class providing the ability to store the factorization of an integer."""

    def __init__(self, **kwargs) -> None:
        if 'value' in kwargs:
            self.__factors = factor(kwargs['value'])
        if 'factors' in kwargs:
            factors = kwargs['factors']
            self.__validate(factors)
            self.__factors = self.__to_standard_form(factors)

    def factors(self) -> Dict[int, int]:
        """Function returning the factors of an integer
        in the form of a map from primes to exponents."""
        return self.__factors

    def primes(self) -> List[int]:
        """Function returning the prime factors of an integer."""
        return list(self.__factors.keys())

    def exponents(self) -> List[int]:
        """Function returning the exponents foreach prime factor of an integer."""
        return list(self.__factors.values())

    def has_prime_divisor(self, n: int) -> bool:
        """Returns True if the parameter is a prime factor of the integer."""
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
        """Returns a list of divisors of an integer"""
        return [divisor for divisor in DivisorsIterator(self)]

    def divisors_count(self) -> int:
        """Returns the number of divisors"""
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


class DivisorsIterator:
    """Class providing the ability to iterate through the divisors of a factored integer."""

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
