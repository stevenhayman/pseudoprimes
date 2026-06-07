#include "pseudoprime/large_carmichael/weight.hpp"

namespace pseudoprime {

std::size_t weight(const Integer& a, const std::vector<Integer>& components) {
  for (std::size_t i = components.size(); i-- > 0;) {
    if (a % components[i] != 1) return i + 1;  // 1-based index of the top nontrivial component
  }
  return 0;
}

}  // namespace pseudoprime
