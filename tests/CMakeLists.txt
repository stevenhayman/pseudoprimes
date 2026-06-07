# Prefer a GoogleTest provided by the toolchain (vcpkg, system, MSYS2/pacman);
# fall back to fetching it if absent.
find_package(GTest CONFIG QUIET)
if(NOT GTest_FOUND)
  find_package(GTest MODULE QUIET)
endif()

if(NOT GTest_FOUND)
  message(STATUS "GoogleTest not found locally; fetching it via FetchContent.")
  include(FetchContent)
  FetchContent_Declare(googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

add_executable(pseudoprime_tests
  factorization_test.cpp
  modular_integer_test.cpp
  modular_product_test.cpp
  subset_product_test.cpp
  primality_test.cpp
  carmichael_test.cpp
  weight_test.cpp
  subgroup_chain_test.cpp
  large_carmichael_test.cpp
  layered_solver_test.cpp
  erdos_construction_test.cpp
  sieve_prime_set_test.cpp
  highly_composite_test.cpp)
target_link_libraries(pseudoprime_tests PRIVATE pseudoprime GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(pseudoprime_tests)
