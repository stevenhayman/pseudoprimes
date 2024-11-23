"""Module for integer-related mathematical operations.

This module provides fundamental mathematical operations for working with integers,
including logarithms, modular arithmetic, and integer roots.
"""

def int_log(b: int, n: int) -> int:
    """Calculate the integer logarithm (floor) of n base b.
    
    Args:
        b: The base of the logarithm
        n: The number to calculate the logarithm of
        
    Returns:
        int: The floor of log base b of n
        
    Example:
        >>> int_log(2, 7)
        2    # since 2^2 = 4 < 7 < 8 = 2^3
    """
    i = 0
    value = 1
    while value <= n:
        value *= b
        i += 1
    return i - 1

def sub_mod(a: int, b: int, n: int) -> int:
    """Compute the modular subtraction of two numbers.
    
    Args:
        a: The minuend
        b: The subtrahend
        n: The modulus
        
    Returns:
        int: (a - b) mod n, guaranteed to be in the range [0, n)
        
    Example:
        >>> sub_mod(2, 7, 5)
        0    # because (2 - 7) ≡ 0 (mod 5)
    """
    t = (a % n) - (b % n)
    if t < 0:
        t += n
    elif t > n:
        t -= n
    return t

def int_root(n: int, k: int) -> int:
    """Calculate the integer k-th root (floor) of n.
    
    Uses Newton's method to compute the floor of the k-th root of n.
    
    Args:
        n: The number to find the root of
        k: The root degree (e.g., 2 for square root, 3 for cube root)
        
    Returns:
        int: The floor of the k-th root of n
        
    Example:
        >>> int_root(10, 2)
        3    # since 3^2 = 9 < 10 < 16 = 4^2
    """
    u = n
    while True:
        s = u
        t = (k - 1) * s + n // (s ** (k - 1))
        u = t // k
        if u >= s:
            break
    return s
