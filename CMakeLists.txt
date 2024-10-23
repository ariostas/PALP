cmake_minimum_required(VERSION 3.14...3.30)

project(PALP
  VERSION
    2.21
  DESCRIPTION
    "A Package for Analyzing Lattice Polytopes "
  LANGUAGES
    C
)

set(POLY_Dmax 6 CACHE STRING "Maximum polytope dimension")

# Compiling with -DNDEBUG breaks the code
string(REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCES Coord.c Rat.c Vertex.c Polynf.c)
set(CLASS_SRC Subpoly.c Subadd.c Subdb.c)
set(NEF_SRC E_Poly.c Nefpart.c LG.c)
set(MORI_SRC MoriCone.c SingularInput.c LG.c)

add_library(objects STATIC ${SOURCES})
target_compile_definitions(objects PUBLIC POLY_Dmax=${POLY_Dmax})
add_library(class_objects STATIC ${CLASS_SRC})
target_compile_definitions(class_objects PUBLIC POLY_Dmax=${POLY_Dmax})
add_library(nef_objects STATIC ${NEF_SRC})
target_compile_definitions(nef_objects PUBLIC POLY_Dmax=${POLY_Dmax})
add_library(mori_objects STATIC ${MORI_SRC})
target_compile_definitions(mori_objects PUBLIC POLY_Dmax=${POLY_Dmax})

include_directories(${CMAKE_SOURCE_DIR})

add_executable(poly.x poly.c LG.c)
target_link_libraries(poly.x PRIVATE objects)
target_compile_definitions(poly.x PRIVATE POLY_Dmax=${POLY_Dmax})

add_executable(class.x class.c)
target_link_libraries(class.x PRIVATE objects class_objects)
target_compile_definitions(class.x PRIVATE POLY_Dmax=${POLY_Dmax})

add_executable(cws.x cws.c LG.c)
target_link_libraries(cws.x PRIVATE objects)
target_compile_definitions(cws.x PRIVATE POLY_Dmax=${POLY_Dmax})

add_executable(nef.x nef.c)
target_link_libraries(nef.x PRIVATE objects nef_objects)
target_compile_definitions(nef.x PRIVATE POLY_Dmax=${POLY_Dmax})

add_executable(mori.x mori.c)
target_link_libraries(mori.x PRIVATE objects mori_objects)
target_compile_definitions(mori.x PRIVATE POLY_Dmax=${POLY_Dmax})
