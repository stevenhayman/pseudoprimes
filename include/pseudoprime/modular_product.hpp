#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "pseudoprime/modular_integer.hpp"

namespace pseudoprime {

// A product of chosen elements reduced modulo a fixed modulus, carried together
// with the indices of the elements that formed it. The residue and its modulus
// are held by the composed ModularInteger, so no separate modulus field is
// needed.
//
// This is the recurring "a subset's product, plus which elements" concept shared
// by combining search strategies (the large-Carmichael pairing engine,
// meet-in-the-middle). It is deliberately not used by depth-first backtracking
// search, where mutating one accumulator with push/pop is cheaper than copying
// the element list per branch.
class ModularProduct {
 public:
  // The empty product (residue 1) modulo `modulus`.
  explicit ModularProduct(Integer modulus) : value_(std::move(modulus)) {}

  ModularProduct(ModularInteger value, std::vector<std::size_t> elements)
      : value_(std::move(value)), elements_(std::move(elements)) {}

  const Integer& residue() const { return value_.value(); }
  const Integer& modulus() const { return value_.modulus(); }
  const std::vector<std::size_t>& elements() const { return elements_; }
  std::size_t size() const { return elements_.size(); }

  // Multiply in a single element (its residue) recorded under `index`.
  void multiply(const Integer& element_residue, std::size_t index) {
    value_ *= element_residue;
    elements_.push_back(index);
  }

  // Combine with another product over the same modulus: residues multiply,
  // element lists concatenate (this product's elements first, then `other`'s).
  ModularProduct combine(const ModularProduct& other) const {
    ModularInteger v = value_;
    v *= other.value_;
    std::vector<std::size_t> e = elements_;
    e.insert(e.end(), other.elements_.begin(), other.elements_.end());
    return ModularProduct(std::move(v), std::move(e));
  }

 private:
  ModularInteger value_;
  std::vector<std::size_t> elements_;
};

}  // namespace pseudoprime
