# FindGMP.cmake
#
# Locates the GNU Multiple Precision Arithmetic Library (GMP) and its C++
# interface (gmpxx). GMP ships no CMake config file, so we locate the headers
# and libraries ourselves, following the modern CMake find-module idiom:
# pkg-config supplies hints where available (MSYS2/pacman and most Linux
# installs ship gmp.pc / gmpxx.pc), and find_path/find_library confirm the
# actual paths and degrade gracefully where it does not (some vcpkg layouts).
# The GMP_ROOT variable / environment variable is honoured for manual installs.
#
# Provides imported targets:
#   GMP::gmp    - the C library
#   GMP::gmpxx  - the C++ wrapper (mpz_class etc.); links GMP::gmp
#
# Sets:
#   GMP_FOUND, GMP_VERSION,
#   GMP_INCLUDE_DIR, GMPXX_INCLUDE_DIR, GMP_LIBRARY, GMPXX_LIBRARY

# pkg-config is only a hint provider here; its absence is not fatal.
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(PC_GMP   QUIET gmp)
  pkg_check_modules(PC_GMPXX QUIET gmpxx)
endif()

find_path(GMP_INCLUDE_DIR
  NAMES gmp.h
  HINTS ${PC_GMP_INCLUDE_DIRS} ${GMP_ROOT} ENV GMP_ROOT
  PATH_SUFFIXES include
  DOC "GMP include directory")

find_path(GMPXX_INCLUDE_DIR
  NAMES gmpxx.h
  HINTS ${PC_GMPXX_INCLUDE_DIRS} ${GMP_ROOT} ENV GMP_ROOT
  PATH_SUFFIXES include
  DOC "GMP C++ include directory")

find_library(GMP_LIBRARY
  NAMES gmp libgmp
  HINTS ${PC_GMP_LIBRARY_DIRS} ${GMP_ROOT} ENV GMP_ROOT
  PATH_SUFFIXES lib
  DOC "GMP C library")

find_library(GMPXX_LIBRARY
  NAMES gmpxx libgmpxx
  HINTS ${PC_GMPXX_LIBRARY_DIRS} ${GMP_ROOT} ENV GMP_ROOT
  PATH_SUFFIXES lib
  DOC "GMP C++ library")

# Derive the version from the public macros in gmp.h.
if(GMP_INCLUDE_DIR AND EXISTS "${GMP_INCLUDE_DIR}/gmp.h")
  file(STRINGS "${GMP_INCLUDE_DIR}/gmp.h" _gmp_version_lines
    REGEX "^#define __GNU_MP_VERSION(_MINOR|_PATCHLEVEL)?[ \t]+[0-9]+")
  foreach(_line IN LISTS _gmp_version_lines)
    if(_line MATCHES "__GNU_MP_VERSION[ \t]+([0-9]+)")
      set(_gmp_major "${CMAKE_MATCH_1}")
    elseif(_line MATCHES "__GNU_MP_VERSION_MINOR[ \t]+([0-9]+)")
      set(_gmp_minor "${CMAKE_MATCH_1}")
    elseif(_line MATCHES "__GNU_MP_VERSION_PATCHLEVEL[ \t]+([0-9]+)")
      set(_gmp_patch "${CMAKE_MATCH_1}")
    endif()
  endforeach()
  if(DEFINED _gmp_major AND DEFINED _gmp_minor AND DEFINED _gmp_patch)
    set(GMP_VERSION "${_gmp_major}.${_gmp_minor}.${_gmp_patch}")
  endif()
  unset(_gmp_version_lines)
  unset(_gmp_major)
  unset(_gmp_minor)
  unset(_gmp_patch)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
  REQUIRED_VARS GMP_LIBRARY GMPXX_LIBRARY GMP_INCLUDE_DIR GMPXX_INCLUDE_DIR
  VERSION_VAR GMP_VERSION)

if(GMP_FOUND)
  if(NOT TARGET GMP::gmp)
    add_library(GMP::gmp UNKNOWN IMPORTED)
    set_target_properties(GMP::gmp PROPERTIES
      IMPORTED_LOCATION "${GMP_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}")
  endif()

  if(NOT TARGET GMP::gmpxx)
    add_library(GMP::gmpxx UNKNOWN IMPORTED)
    set_target_properties(GMP::gmpxx PROPERTIES
      IMPORTED_LOCATION "${GMPXX_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${GMPXX_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES GMP::gmp)
  endif()
endif()

mark_as_advanced(GMP_INCLUDE_DIR GMPXX_INCLUDE_DIR GMP_LIBRARY GMPXX_LIBRARY)
