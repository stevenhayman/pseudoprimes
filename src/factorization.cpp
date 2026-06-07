#include "pseudoprime/factorization.hpp"

#include <algorithm>
#include <generator>

namespace pseudoprime {

Factorization::Factorization(std::vector<PrimePower> factors)
    : factors_(std::move(factors)) {
  std::sort(factors_.begin(), factors_.end(),
            [](const PrimePower& a, const PrimePower& b) { return a.prime < b.prime; });
}

Factorization Factorization::from_prime_power(const Integer& prime, unsigned long exponent) {
  if (exponent == 0) return Factorization{};
  return Factorization(std::vector<PrimePower>{{prime, exponent}});
}

Factorization Factorization::from_primes(const std::vector<Integer>& primes) {
  std::vector<Integer> sorted(primes);
  std::sort(sorted.begin(), sorted.end());

  std::vector<PrimePower> factors;
  for (const Integer& p : sorted) {
    if (!factors.empty() && factors.back().prime == p) {
      ++factors.back().exponent;
    } else {
      factors.push_back({p, 1});
    }
  }
  return Factorization(std::move(factors));
}

Integer Factorization::value() const {
  Integer n = 1;
  for (const PrimePower& pe : factors_) n *= pow_ui(pe.prime, pe.exponent);
  return n;
}

Integer Factorization::num_divisors() const {
  Integer d = 1;
  for (const PrimePower& pe : factors_) d *= (pe.exponent + 1);
  return d;
}

bool Factorization::is_squarefree() const {
  for (const PrimePower& pe : factors_) {
    if (pe.exponent != 1) return false;
  }
  return true;
}

std::vector<Integer> Factorization::primes() const {
  std::vector<Integer> out;
  out.reserve(factors_.size());
  for (const PrimePower& pe : factors_) out.push_back(pe.prime);
  return out;
}

std::vector<Integer> Factorization::prime_powers() const {
  std::vector<Integer> out;
  out.reserve(factors_.size());
  for (const PrimePower& pe : factors_) out.push_back(pow_ui(pe.prime, pe.exponent));
  return out;  // factors_ is sorted ascending by prime
}

unsigned long Factorization::multiplicity(const Integer& p) const {
  for (const PrimePower& pe : factors_) {
    if (pe.prime == p) return pe.exponent;
  }
  return 0;
}

Integer Factorization::euler_phi() const {
  // phi(n) = product over p^e of p^(e-1) * (p - 1).
  Integer phi = 1;
  for (const PrimePower& pe : factors_) {
    phi *= pow_ui(pe.prime, pe.exponent - 1) * (pe.prime - 1);
  }
  return phi;
}

Integer Factorization::carmichael_lambda() const {
  // lambda(n) = lcm over prime powers of lambda(p^e), where
  //   lambda(p^e) = p^(e-1) * (p - 1)   for odd p (or p=2, e<=2),
  //   lambda(2^e) = 2^(e-2)             for e >= 3.
  Integer lambda = 1;
  for (const PrimePower& pe : factors_) {
    Integer term;
    if (pe.prime == 2 && pe.exponent >= 3) {
      term = pow_ui(Integer(2), pe.exponent - 2);
    } else {
      term = pow_ui(pe.prime, pe.exponent - 1) * (pe.prime - 1);
    }
    lambda = lcm(lambda, term);
  }
  return lambda;
}

std::vector<Integer> Factorization::divisors() const {
  std::vector<Integer> divs{Integer(1)};
  for (const PrimePower& pe : factors_) {
    const std::size_t base_count = divs.size();
    Integer power = 1;
    for (unsigned long e = 1; e <= pe.exponent; ++e) {
      power *= pe.prime;
      for (std::size_t i = 0; i < base_count; ++i) divs.push_back(divs[i] * power);
    }
  }
  std::sort(divs.begin(), divs.end());
  return divs;
}

namespace {

// DFS over the per-prime exponent choices (0..e), holding only the divisor
// currently being built; yields it at each leaf. `current` lives in the caller's
// coroutine frame and is shared by reference across the recursion.
std::generator<Factorization> generate_divisors(const std::vector<PrimePower>& factors,
                                                std::size_t i, std::vector<PrimePower>& current) {
  if (i == factors.size()) {
    co_yield Factorization(current);  // current is built in ascending-prime order
    co_return;
  }
  co_yield std::ranges::elements_of(generate_divisors(factors, i + 1, current));  // omit prime i
  const PrimePower& pe = factors[i];
  for (unsigned long e = 1; e <= pe.exponent; ++e) {
    current.push_back({pe.prime, e});
    co_yield std::ranges::elements_of(generate_divisors(factors, i + 1, current));
    current.pop_back();
  }
}

}  // namespace

std::generator<Factorization> Factorization::divisors_lazy() const {
  std::vector<PrimePower> current;
  co_yield std::ranges::elements_of(generate_divisors(factors_, 0, current));
}

std::vector<Factorization> Factorization::divisors_factored() const {
  std::vector<Factorization> out;
  for (const Factorization& d : divisors_lazy()) out.push_back(d);
  std::sort(out.begin(), out.end(), [](const Factorization& a, const Factorization& b) {
    return a.value() < b.value();
  });
  return out;
}

Factorization gcd(const Factorization& a, const Factorization& b) {
  const std::vector<PrimePower>& fa = a.factors();
  const std::vector<PrimePower>& fb = b.factors();

  std::vector<PrimePower> result;
  std::size_t i = 0, j = 0;
  while (i < fa.size() && j < fb.size()) {
    if (fa[i].prime < fb[j].prime) {
      ++i;
    } else if (fb[j].prime < fa[i].prime) {
      ++j;
    } else {
      result.push_back({fa[i].prime, std::min(fa[i].exponent, fb[j].exponent)});
      ++i;
      ++j;
    }
  }
  return Factorization(std::move(result));
}

}  // namespace pseudoprime
