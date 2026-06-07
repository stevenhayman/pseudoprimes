#pragma once

#include "pseudoprime/integer.hpp"

namespace pseudoprime {

// An integer considered modulo a fixed modulus m: an element of Z/mZ.
//
// The value is always kept as the canonical representative in [0, m); the
// original (pre-reduction) operand is never retained. So this is effectively a
// residue-class representative paired with its modulus -- "ModularInteger"
// names the element we manipulate, while the modulus identifies which ring
// Z/mZ it lives in. Centralizing reduction here means callers (the
// subset-product accumulator in particular) cannot forget to reduce.
//
// The default-constructed value is the multiplicative identity 1. All operands
// of a binary operation must share the same modulus.
class ModularInteger {
 public:
  explicit ModularInteger(Integer modulus) : mod_(std::move(modulus)), val_(1) {}

  ModularInteger(Integer value, Integer modulus)
      : mod_(std::move(modulus)), val_(std::move(value)) {
    normalize();
  }

  const Integer& value() const { return val_; }
  const Integer& modulus() const { return mod_; }

  ModularInteger& operator*=(const Integer& rhs) {
    val_ *= (rhs % mod_);
    normalize();
    return *this;
  }

  ModularInteger& operator*=(const ModularInteger& rhs) {
    val_ *= rhs.val_;
    normalize();
    return *this;
  }

  friend ModularInteger operator*(ModularInteger a, const ModularInteger& b) {
    a *= b;
    return a;
  }

  bool operator==(const Integer& rhs) const {
    Integer r = rhs % mod_;
    if (r < 0) r += mod_;
    return val_ == r;
  }

 private:
  void normalize() {
    val_ %= mod_;
    if (val_ < 0) val_ += mod_;
  }

  Integer mod_;
  Integer val_;
};

}  // namespace pseudoprime
