cmake_minimum_required(VERSION 3.14...4.1)

project(PALP
  VERSION
    2.21
  DESCRIPTION
    "A Package for Analyzing Lattice Polytopes "
  LANGUAGES
    C
)

# PALP is built both standalone and as a subproject (the Python bindings pull
# it in with add_subdirectory()). Executables, install rules and tests are only
# of interest to the standalone build.
if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(PALP_IS_TOP_LEVEL TRUE)
else()
  set(PALP_IS_TOP_LEVEL FALSE)
endif()

set(POLY_Dmax 6 CACHE STRING "Maximum polytope dimension")
if (NOT "${POLY_Dmax}" MATCHES "^[1-9][0-9]*$")
  message(FATAL_ERROR "POLY_Dmax must be a positive integer, got '${POLY_Dmax}'")
endif()

get_property(PALP_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

# The GNUmakefile compiles with -O3 by default. Without this, a plain
# "cmake -S . -B build" would produce unoptimized binaries, which for PALP
# is a difference of several orders of magnitude in run time.
if (PALP_IS_TOP_LEVEL AND NOT PALP_MULTI_CONFIG AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release RelWithDebInfo MinSizeRel)
endif()

# PALP calls functions with side effects from inside assert(), e.g.
#   assert( 1 == W_to_GLZ(WZ->w,&a,B) );   [Polynf.c]
# so defining NDEBUG silently removes real work. CMake adds -DNDEBUG to every
# optimized configuration, so strip it from all of them.
foreach(config RELEASE RELWITHDEBINFO MINSIZEREL)
  string(REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
  string(REPLACE "/DNDEBUG" "" CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
endforeach()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# CPPFLAGS from the GNUmakefile. _FILE_OFFSET_BITS=64 changes the width of
# off_t and _LARGEFILE_SOURCE additionally switches Subpoly.h from
# fseek()/ftell() to fseeko()/ftello(), so both have to be seen by every
# translation unit that is linked together -- hence PUBLIC below.
set(PALP_COMPILE_DEFINITIONS _FILE_OFFSET_BITS=64 _LARGEFILE_SOURCE)

set(SOURCES Coord.c Rat.c Vertex.c Polynf.c)
set(CLASS_SRC Subpoly.c Subadd.c Subdb.c)
set(NEF_SRC E_Poly.c Nefpart.c LG.c)
set(MORI_SRC MoriCone.c SingularInput.c LG.c)

# Every object file has to be compiled with the POLY_Dmax it will be linked
# with: Global.h derives array sizes and struct layouts from it, so mixing
# values within one executable is undefined behaviour. That is why the
# libraries below are created once per dimension rather than shared.
function(palp_add_library target dmax)
  # Sources are named relative to the project root so that these helpers also
  # work when called from a subdirectory (see tests/CMakeLists.txt).
  set(sources "")
  foreach(source IN LISTS ARGN)
    list(APPEND sources "${PROJECT_SOURCE_DIR}/${source}")
  endforeach()
  add_library(${target} STATIC ${sources})
  target_compile_definitions(${target}
    PUBLIC POLY_Dmax=${dmax} ${PALP_COMPILE_DEFINITIONS})
  target_include_directories(${target} PUBLIC ${PROJECT_SOURCE_DIR})
endfunction()

# Creates objects${suffix}, class_objects${suffix}, nef_objects${suffix} and
# mori_objects${suffix}, all compiled with POLY_Dmax=${dmax}.
function(palp_add_libraries suffix dmax)
  palp_add_library("objects${suffix}"       ${dmax} ${SOURCES})
  palp_add_library("class_objects${suffix}" ${dmax} ${CLASS_SRC})
  palp_add_library("nef_objects${suffix}"   ${dmax} ${NEF_SRC})
  palp_add_library("mori_objects${suffix}"  ${dmax} ${MORI_SRC})
endfunction()

# Creates poly${suffix}.x, class${suffix}.x, cws${suffix}.x, nef${suffix}.x and
# mori${suffix}.x from the libraries created by palp_add_libraries(${suffix}).
# The resulting target names are returned in PALP_PROGRAM_TARGETS.
function(palp_add_programs suffix)
  add_executable(poly${suffix}.x ${PROJECT_SOURCE_DIR}/poly.c ${PROJECT_SOURCE_DIR}/LG.c)
  target_link_libraries(poly${suffix}.x PRIVATE objects${suffix})

  add_executable(class${suffix}.x ${PROJECT_SOURCE_DIR}/class.c)
  target_link_libraries(class${suffix}.x PRIVATE objects${suffix} class_objects${suffix})

  add_executable(cws${suffix}.x ${PROJECT_SOURCE_DIR}/cws.c ${PROJECT_SOURCE_DIR}/LG.c)
  target_link_libraries(cws${suffix}.x PRIVATE objects${suffix})

  add_executable(nef${suffix}.x ${PROJECT_SOURCE_DIR}/nef.c)
  target_link_libraries(nef${suffix}.x PRIVATE objects${suffix} nef_objects${suffix})

  add_executable(mori${suffix}.x ${PROJECT_SOURCE_DIR}/mori.c)
  target_link_libraries(mori${suffix}.x PRIVATE objects${suffix} mori_objects${suffix})

  set(PALP_PROGRAM_TARGETS
    poly${suffix}.x class${suffix}.x cws${suffix}.x nef${suffix}.x mori${suffix}.x
    PARENT_SCOPE)
endfunction()

palp_add_libraries("" ${POLY_Dmax})

if (PALP_IS_TOP_LEVEL)
  palp_add_programs("")

  include(GNUInstallDirs)
  install(TARGETS ${PALP_PROGRAM_TARGETS}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  option(PALP_BUILD_TESTS
    "Build the PALP test suite (compiles one set of programs per test dimension)"
    OFF)
  if (PALP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
  endif()
endif()
