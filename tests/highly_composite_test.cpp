#include "pseudoprime/highly_composite.hpp"

#include <gtest/gtest.h>

#include <set>
#include <vector>

#include "pseudoprime/erdos_construction.hpp"
#include "pseudoprime/factorization.hpp"

namespace {

using pseudoprime::Factorization;
using pseudoprime::Integer;

// Collect the first `n` terms of a Factorization generator.
template <typename Range>
std::vector<Factorization> take(Range&& range, std::size_t n) {
  std::vector<Factorization> out;
  for (const Factorization& f : range) {
    if (out.size() >= n) break;
    out.push_back(f);
  }
  return out;
}

std::vector<Integer> values(const std::vector<Factorization>& fs) {
  std::vector<Integer> out;
  out.reserve(fs.size());
  for (const Factorization& f : fs) out.push_back(f.value());
  return out;
}

TEST(FactorizationFromExponents, BuildsAndRoundTrips) {
  // 2^3 * 3^2 * 5 * 7 = 2520
  const Factorization f = pseudoprime::factorization_from_exponents({3, 2, 1, 1});
  EXPECT_EQ(f.value(), Integer(2520));
  ASSERT_EQ(f.num_distinct_primes(), 4u);
  EXPECT_EQ(f.multiplicity(Integer(2)), 3u);
  EXPECT_EQ(f.multiplicity(Integer(3)), 2u);
}

TEST(FactorizationFromExponents, DropsZeroExponents) {
  // 2^2 * 5 = 20 (the 3 slot is zero)
  const Factorization f = pseudoprime::factorization_from_exponents({2, 0, 1});
  EXPECT_EQ(f.value(), Integer(20));
  EXPECT_EQ(f.num_distinct_primes(), 2u);
  EXPECT_EQ(f.multiplicity(Integer(3)), 0u);
}

TEST(NumDivisors, MatchesDivisorEnumeration) {
  const std::vector<Factorization> cases = {
      Factorization{},  // 1 -> d = 1
      Factorization({{Integer(2), 4}, {Integer(5), 1}}),               // 80 -> 10
      pseudoprime::factorization_from_exponents({3, 2, 1, 1}),         // 2520 -> 48
      pseudoprime::factorization_from_exponents({2, 1, 1}),            // 60 -> 12
  };
  for (const Factorization& f : cases) {
    EXPECT_EQ(f.num_divisors(), Integer(static_cast<unsigned long>(f.divisors().size())))
        << "n = " << f.value();
  }
}

TEST(SuperiorHighlyComposite, MatchesA002201) {
  const std::vector<Integer> expected = {2,    6,    12,    60,     120,
                                         360,  2520, 5040,  55440,  720720};
  EXPECT_EQ(values(take(pseudoprime::superior_highly_composite(), expected.size())), expected);
}

TEST(SuperiorHighlyComposite, CountBounds) {
  EXPECT_EQ(take(pseudoprime::superior_highly_composite(4), 100).size(), 4u);
}

TEST(SuperiorHighlyComposite, FactorizationOf2520) {
  // 7th term is 2520 = 2^3 * 3^2 * 5 * 7.
  const std::vector<Factorization> fs = take(pseudoprime::superior_highly_composite(7), 7);
  const Factorization& f = fs.back();
  EXPECT_EQ(f.value(), Integer(2520));
  EXPECT_EQ(f.multiplicity(Integer(2)), 3u);
  EXPECT_EQ(f.multiplicity(Integer(3)), 2u);
  EXPECT_EQ(f.multiplicity(Integer(5)), 1u);
  EXPECT_EQ(f.multiplicity(Integer(7)), 1u);
}

TEST(CandidateHighlyComposite, ValueOrderedAndNonIncreasing) {
  const std::vector<Factorization> fs = take(pseudoprime::candidate_highly_composite(), 200);
  std::set<Integer> distinct;
  for (std::size_t i = 0; i < fs.size(); ++i) {
    // Strictly increasing value (also proves no duplicates).
    if (i > 0) EXPECT_LT(fs[i - 1].value(), fs[i].value());
    distinct.insert(fs[i].value());

    // Exponents non-increasing across the prime run.
    const auto& factors = fs[i].factors();
    for (std::size_t j = 1; j < factors.size(); ++j) {
      EXPECT_GE(factors[j - 1].exponent, factors[j].exponent) << "n = " << fs[i].value();
    }
  }
  EXPECT_EQ(distinct.size(), fs.size());
}

TEST(HighlyComposite, MatchesA002182) {
  const std::vector<Integer> expected = {1,   2,   4,   6,   12,  24,  36,   48,   60,    120,
                                         180, 240, 360, 720, 840, 1260, 1680, 2520, 5040};
  EXPECT_EQ(values(take(pseudoprime::highly_composite(), expected.size())), expected);
}

TEST(HighlyComposite, DivisorCountStrictlyIncreasing) {
  const std::vector<Factorization> fs = take(pseudoprime::highly_composite(), 25);
  for (std::size_t i = 1; i < fs.size(); ++i) {
    EXPECT_LT(fs[i - 1].num_divisors(), fs[i].num_divisors());
  }
}

TEST(HighlyComposite, ContainsSuperiorSubsequence) {
  const std::vector<Factorization> hcn = take(pseudoprime::highly_composite(), 40);
  std::set<Integer> hcn_values;
  for (const Factorization& f : hcn) hcn_values.insert(f.value());

  // Every SHCN value within the HCN range we sampled must appear among the HCNs.
  const Integer ceiling = hcn.back().value();
  for (const Factorization& s : pseudoprime::superior_highly_composite()) {
    if (s.value() > ceiling) break;
    EXPECT_TRUE(hcn_values.count(s.value()) > 0) << "SHCN " << s.value() << " missing from HCN";
  }
}

TEST(HighlyComposite, EveryHcnIsCandidate) {
  std::set<Integer> candidates;
  for (const Factorization& f : take(pseudoprime::candidate_highly_composite(), 300)) {
    candidates.insert(f.value());
  }
  for (const Factorization& h : take(pseudoprime::highly_composite(), 15)) {
    if (h.value() == 1) continue;  // 1 is not in the candidate stream (starts at 2)
    EXPECT_TRUE(candidates.count(h.value()) > 0) << "HCN " << h.value() << " missing from candidates";
  }
}

TEST(HighlyComposite, FeedsPrimeSet) {
  // A mid-size SHCN used as L: P(L) is nonempty and bounded by d(L).
  const std::vector<Factorization> fs = take(pseudoprime::superior_highly_composite(7), 7);
  const Factorization& L = fs.back();  // 2520
  const auto primes = pseudoprime::generate_prime_set(L);
  EXPECT_GT(primes.size(), 0u);
  EXPECT_LE(Integer(static_cast<unsigned long>(primes.size())), L.num_divisors());
}

}  // namespace
