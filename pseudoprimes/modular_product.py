"""Module for handling modular arithmetic products.

This module provides the ModularProduct class for performing modular arithmetic
operations on products of integers. It maintains both the factors and their
product modulo n, allowing efficient multiplication operations while preserving
the factor history.
"""
from __future__ import annotations

from typing import List, Optional


def multiply_modulo(factors: List[int], modulus: int) -> int:
    """Calculate the product of integers modulo a given number.
    
    Args:
        factors: List of integers to multiply together
        modulus: The modulus to apply to the product
        
    Returns:
        int: The product of all factors modulo the given modulus
        
    Example:
        >>> multiply_modulo([2, 3, 4], 5)
        4    # because (2 * 3 * 4) mod 5 = 24 mod 5 = 4
    """
    value = 1
    for factor in factors:
        value = (value * factor) % modulus
    return value


class ModularProduct:
    """A class for handling products in modular arithmetic.
    
    This class maintains both the individual factors and their product modulo n,
    allowing for efficient multiplication while preserving the history of factors
    that contributed to the product.
    
    Attributes:
        modulus: The modulus for all arithmetic operations
        factors: The list of factors in the product
    """

    def __init__(self, modulus: int, factors: Optional[List[int]] = None) -> None:
        """Initialize a ModularProduct instance.
        
        Args:
            modulus: The modulus for all arithmetic operations
            factors: Optional list of initial factors. Defaults to empty list if None
        """
        self.__modulus = modulus
        self.__factors = factors or []
        self.__value = multiply_modulo(factors, modulus)

    def __mul__(self, arg: any) -> ModularProduct:
        """Multiply this ModularProduct with another value.
        
        Args:
            arg: Either an integer or another ModularProduct instance
            
        Returns:
            ModularProduct: A new instance containing the product
            
        Raises:
            NotImplementedError: If arg is neither an int nor ModularProduct
            ValueError: If multiplying with a ModularProduct with different modulus
        """
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
        """Create a deep copy of this ModularProduct instance.
        
        Returns:
            ModularProduct: A new instance with the same modulus and copied factors
        """
        return ModularProduct(self.__modulus, self.__factors.copy())

    def value(self) -> int:
        """Get the current value of the product modulo n.
        
        Returns:
            int: The product of all factors modulo the instance's modulus
        """
        return self.__value

    def factors(self) -> List[int]:
        """Get the list of factors in this product.
        
        Returns:
            List[int]: A list of all factors that have been multiplied together
        """
        return self.__factors

    def modulus(self) -> int:
        """Retrieves the modulus of the product."""
        return self.__modulus
