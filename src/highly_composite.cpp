#include "pseudoprime/highly_composite.hpp"

#include <cmath>
#include <queue>
#include <set>
#include <utility>
#include <vector>

namespace pseudoprime {

namespace {

// Trial-division successor in the prime sequence. The primes needed here are
// tiny (well under 200 for any reachable size), so a deterministic check keeps
// us clear of any probable-prime concern -- in keeping with the library's
// no-probable-primes stance -- at negligible cost.
bool is_prime_small(const Integer& n) {
  if (n < 2) return false;
  for (Integer d = 2; d * d <= n; ++d) {
    if (n % d == 0) return false;
  }
  return true;
}

Integer next_prime(const Integer& p) {
  Integer q = p + 1;
  while (!is_prime_small(q)) ++q;
  return q;
}

// Extend `primes` (2, 3, 5, ...) until it holds at least `count` entries.
void ensure_primes(std::vector<Integer>& primes, std::size_t count) {
  if (primes.empty()) primes.push_back(Integer(2));
  while (primes.size() < count) primes.push_back(next_prime(primes.back()));
}

// Assemble a Factorization from a per-index exponent vector over `primes`,
// dropping zero exponents. `primes` must cover every nonzero index.
Factorization build(const std::vector<Integer>& primes,
                    const std::vector<unsigned long>& exps) {
  std::vector<PrimePower> factors;
  for (std::size_t i = 0; i < exps.size(); ++i) {
    if (exps[i] != 0) factors.push_back({primes[i], exps[i]});
  }
  return Factorization(std::move(factors));
}

}  // namespace

Factorization factorization_from_exponents(const std::vector<unsigned long>& exponents) {
  std::vector<Integer> primes;
  ensure_primes(primes, exponents.size());
  return build(primes, exponents);
}

std::generator<Factorization> superior_highly_composite(std::size_t count) {
  // One event per prime: the eps threshold at which that prime's exponent next
  // increases. crit(p, a) = ln((a+1)/a) / ln(p) is the eps below which the a-th
  // unit of p is present. Sweeping eps downward == popping the largest eps.
  struct Event {
    long double eps;
    std::size_t prime_index;
    unsigned long next_exp;  // the exponent this event raises p to (>= 1)
    bool operator<(const Event& other) const { return eps < other.eps; }  // max-heap
  };

  std::vector<Integer> primes;
  ensure_primes(primes, 1);  // start with prime 2

  auto crit = [&](std::size_t pi, unsigned long a) -> long double {
    const long double ratio = (static_cast<long double>(a) + 1.0L) / static_cast<long double>(a);
    return std::log(ratio) / std::log(primes[pi].get_d());
  };

  std::priority_queue<Event> heap;
  heap.push({crit(0, 1), 0, 1});  // prime 2, first unit; crit == 1

  std::vector<unsigned long> exps;  // exponent per prime index
  std::size_t produced = 0;

  while (count == 0 || produced < count) {
    const Event ev = heap.top();
    heap.pop();

    if (ev.prime_index >= exps.size()) exps.resize(ev.prime_index + 1, 0);
    exps[ev.prime_index] = ev.next_exp;  // raise p's exponent by one

    co_yield build(primes, exps);
    ++produced;

    // Next unit of the same prime.
    heap.push({crit(ev.prime_index, ev.next_exp + 1), ev.prime_index, ev.next_exp + 1});

    // The first use of a prime unlocks the next prime (whose first-unit eps is
    // strictly smaller, so it is seeded well before it is needed).
    if (ev.next_exp == 1) {
      const std::size_t next_index = ev.prime_index + 1;
      ensure_primes(primes, next_index + 1);
      heap.push({crit(next_index, 1), next_index, 1});
    }
  }
}

std::generator<Factorization> candidate_highly_composite() {
  // Best-first enumeration of non-increasing exponent vectors in value order.
  struct Node {
    Integer value;
    std::vector<unsigned long> exps;  // non-increasing, all >= 1
    bool operator>(const Node& other) const { return value > other.value; }
  };

  std::vector<Integer> primes;

  std::priority_queue<Node, std::vector<Node>, std::greater<Node>> heap;  // min by value
  std::set<std::vector<unsigned long>> seen;

  ensure_primes(primes, 1);
  const std::vector<unsigned long> start{1};  // n = 2
  heap.push({primes[0], start});
  seen.insert(start);

  auto offer = [&](const Integer& value, std::vector<unsigned long> exps) {
    if (seen.insert(exps).second) heap.push({value, std::move(exps)});
  };

  while (true) {
    const Node node = heap.top();
    heap.pop();

    co_yield build(primes, node.exps);

    const std::size_t k = node.exps.size();

    // (a) Increment exponent i where it keeps the vector non-increasing, i.e.
    //     i == 0 or exps[i-1] > exps[i].
    for (std::size_t i = 0; i < k; ++i) {
      if (i == 0 || node.exps[i - 1] > node.exps[i]) {
        std::vector<unsigned long> child = node.exps;
        ++child[i];
        offer(node.value * primes[i], std::move(child));
      }
    }

    // (b) Append a new prime with exponent 1 (always non-increasing-valid).
    ensure_primes(primes, k + 1);
    std::vector<unsigned long> appended = node.exps;
    appended.push_back(1);
    offer(node.value * primes[k], std::move(appended));
  }
}

std::generator<Factorization> highly_composite() {
  // n=1 (the empty factorization) is the first HCN; d(1)=1. Then keep only the
  // strict divisor-count records of the value-ordered candidate stream.
  co_yield Factorization{};
  Integer best = 1;  // d(1)
  for (const Factorization& cand : candidate_highly_composite()) {
    const Integer d = cand.num_divisors();
    if (d > best) {
      best = d;
      co_yield cand;
    }
  }
}

}  // namespace pseudoprime
