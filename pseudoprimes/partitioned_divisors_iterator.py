"""Iterator for partitioning divisors of a factorization for parallel processing.

This module provides functionality to split the work of iterating through divisors
across multiple processes. It partitions the divisors based on their prime factors
to ensure an even distribution of work.
"""
from __future__ import annotations
from typing import Iterator, List, Tuple

from pseudoprimes.factorization import DivisorsIterator, Factorization
from pseudoprimes.lists import split_equally


class PartitionedDivisorsIterator:
    """Iterator that yields divisors from a specific partition of a factorization.
    
    This class splits the prime factors of a number into two groups:
    1. Partition
    2. Remaining
    """

    def __init__(
        self,
        factorization: Factorization,
        partition_index: int,
        partitions_count: int
    ) -> None:
        self.__factorization = factorization
        self.__partition_index = partition_index
        self.__partitions_count = partitions_count

        (
            self.__partition_keys_factorization,
            self.__remaining_factorization
        ) = self.__get_partition_keys_and_remaining_factorizations()

        self.__partition_key_factorizations = self.__get_partition_key_factorizations()

    def __get_partition_keys_index(self) -> int:
        exponents = self.__factorization.exponents()

        product = 1
        for index, exponent in enumerate(exponents):
            product *= exponent + 1

            if product >= self.__partitions_count:
                return index

    def __get_partition_keys_and_remaining_factorizations(
        self,
    ) -> Tuple[Factorization, Factorization]:
        partition_keys_index = self.__get_partition_keys_index()

        primes = self.__factorization.primes()
        exponents = self.__factorization.exponents()

        partition_keys_factors = dict(zip(
            primes[0:partition_keys_index + 1],
            exponents[0:partition_keys_index + 1]
        ))
        remaining_factors = dict(zip(
            primes[partition_keys_index + 1:],
            exponents[partition_keys_index + 1:]
        ))
        return (
            Factorization(factors=partition_keys_factors),
            Factorization(factors=remaining_factors)
        )

    def __get_partition_key_factorizations(self) -> List[Factorization]:
        partition_keys_divisors = list(
            DivisorsIterator(self.__partition_keys_factorization)
        )

        return split_equally(
            partition_keys_divisors,
            self.__partitions_count
        )[self.__partition_index]

    def __iter__(self) -> Iterator[Factorization]:
        for partition_key_divisor in self.__partition_key_factorizations:
            for divisor in DivisorsIterator(self.__remaining_factorization):
                factors = partition_key_divisor.factors() | divisor.factors()
                yield Factorization(factors=factors)
