cmake_minimum_required(VERSION 3.21)

project(pseudoprimes
  VERSION 0.1.0
  DESCRIPTION "A library for exploring the mathematics of pseudoprime numbers"
  LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(PP_BUILD_TESTS "Build unit tests (requires GoogleTest)" ON)
option(PP_BUILD_APPS  "Build command-line tools"               ON)

# ---------------------------------------------------------------------------
# Global settings
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
# GMP is our single justified third-party dependency (arbitrary precision).
# It is acquired from the platform toolchain (MSYS2/pacman, vcpkg, or system
# package manager) rather than vendored.
find_package(GMP REQUIRED)

# ---------------------------------------------------------------------------
# Core library
# ---------------------------------------------------------------------------
add_library(pseudoprime
  src/factorization.cpp
  src/primality/bls.cpp
  src/primality/fellows_koblitz.cpp
  src/primality/lucas.cpp
  src/primality/primality_test.cpp
  src/order.cpp
  src/subset_product/brute_force.cpp
  src/subset_product/meet_in_middle.cpp
  src/large_carmichael/weight.cpp
  src/large_carmichael/subgroup_chain.cpp
  src/large_carmichael/pairing.cpp
  src/large_carmichael/construction.cpp
  src/large_carmichael/layered_solver.cpp
  src/prime_set_strategy/divisor.cpp
  src/prime_set_strategy/sieve.cpp
  src/erdos_construction.cpp
  src/carmichael.cpp
  src/highly_composite.cpp)

target_include_directories(pseudoprime PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(pseudoprime PUBLIC GMP::gmpxx GMP::gmp)

target_compile_options(pseudoprime PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
  $<$<CXX_COMPILER_ID:MSVC>:/W4>)

# ---------------------------------------------------------------------------
# Apps
# ---------------------------------------------------------------------------
# The apps/ directory holds a local CLI driver that is deliberately NOT version
# controlled (see .gitignore). Build it only when present, so a fresh checkout
# configures cleanly without it.
if(PP_BUILD_APPS AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/apps/CMakeLists.txt")
  add_subdirectory(apps)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(PP_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
