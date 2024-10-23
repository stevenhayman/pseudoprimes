"""Module providing the `ModularProduct` class,
which represents a product of integers under a specified modulus. 
The class supports modular multiplication of integers 
or other `ModularProduct` instances, and keeps track of the 
factors contributing to the product."""
from __future__ import annotations

from typing import List, Optional


def multiply_modulo(factors: List[int], modulus: int,) -> int:
    """Calculates the product of the factors modulo the specified modulus."""
    value = 1
    for factor in factors:
        value = (value * factor) % modulus
    return value


class ModularProduct:
    """A class representing a product of integers modulo a given modulus."""

    def __init__(self, modulus: int, factors: Optional[List[int]] = None) -> None:
        """Initializes a ModularProduct instance with a modulus and an optional list of factors."""
        self.__modulus = modulus
        self.__factors = factors or []
        self.__value = multiply_modulo(factors, modulus)

    def __mul__(self, arg: any) -> ModularProduct:
        """Multiplies the current ModularProduct instance with an integer
        or another ModularProduct instance."""
        if isinstance(arg, int):
            return self.__mul__int(arg)
        elif isinstance(arg, ModularProduct):
            return self.__mul__modular_product(arg)
        else:
            raise NotImplementedError()

    def __mul__int(self, n: int) -> ModularProduct:
        """Multiplies the current product by an integer, returning a new ModularProduct."""
        factors = self.__factors + [n]
        return ModularProduct(self.__modulus, factors)

    def __mul__modular_product(self, other: ModularProduct) -> ModularProduct:
        """Multiplies the current product with another ModularProduct,
        returning a new ModularProduct."""
        if self.modulus() != other.modulus():
            raise ValueError("Moduli must be the same")

        factors = self.factors() + other.factors()
        return ModularProduct(self.modulus(), factors)

    def copy(self) -> ModularProduct:
        """Creates a copy of the current ModularProduct instance."""
        return ModularProduct(self.__modulus, self.__factors.copy())

    def value(self) -> int:
        """Retrieves the value of the modular product."""
        return self.__value

    def factors(self) -> List[int]:
        """Retrieves the list of factors contributing to the product."""
        return self.__factors

    def modulus(self) -> int:
        """Retrieves the modulus of the product."""
        return self.__modulus
