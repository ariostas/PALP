cmake_minimum_required(VERSION 3.14...4.1)

project(PALP
  VERSION
    2.21
  DESCRIPTION
    "A Package for Analyzing Lattice Polytopes "
  LANGUAGES
    CXX
)

set(POLY_Dmax 6 CACHE STRING "Maximum polytope dimension")
set(PALP_TEST_DIMENSIONS "4;5;6;11" CACHE STRING "POLY_Dmax values to build and test with CTest")
set(PALP_TESTS_LONG OFF CACHE BOOL "Run PALP tests that require LONG=1")

# Compiling with -DNDEBUG breaks the code
string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable strict warnings for all C++ targets.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# Public include directories: src/ (legacy private headers) and include/
# (namespaced public headers under palp/).
set(PALP_INCLUDE_DIRS ${PALP_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(PALP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Collect the .cpp version of each PALP source under src/.
function(palp_sources out_var)
  set(result)
  foreach(src ${ARGN})
    set(cpp_src "${PALP_SOURCE_DIR}/${src}.cpp")
    if(EXISTS "${cpp_src}")
      list(APPEND result "${cpp_src}")
    else()
      message(FATAL_ERROR "PALP C++ source not found: ${cpp_src}")
    endif()
  endforeach()
  set(${out_var} ${result} PARENT_SCOPE)
endfunction()

palp_sources(SOURCES Coord Rat Vertex Polynf)
palp_sources(CLASS_SRC Subpoly Subadd Subdb)
palp_sources(NEF_SRC E_Poly Nefpart LG)
palp_sources(MORI_SRC MoriCone SingularInput LG)

add_library(objects STATIC ${SOURCES})
target_compile_definitions(objects PUBLIC POLY_Dmax=${POLY_Dmax})
target_include_directories(objects PUBLIC ${PALP_INCLUDE_DIRS})

add_library(class_objects STATIC ${CLASS_SRC})
target_compile_definitions(class_objects PUBLIC POLY_Dmax=${POLY_Dmax})
target_include_directories(class_objects PUBLIC ${PALP_INCLUDE_DIRS})

add_library(nef_objects STATIC ${NEF_SRC})
target_compile_definitions(nef_objects PUBLIC POLY_Dmax=${POLY_Dmax})
target_include_directories(nef_objects PUBLIC ${PALP_INCLUDE_DIRS})

add_library(mori_objects STATIC ${MORI_SRC})
target_compile_definitions(mori_objects PUBLIC POLY_Dmax=${POLY_Dmax})
target_include_directories(mori_objects PUBLIC ${PALP_INCLUDE_DIRS})

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    include(CTest)

    palp_sources(POLY_SRC poly LG)
    add_executable(poly.x ${POLY_SRC})
    target_link_libraries(poly.x PRIVATE objects)
    target_compile_definitions(poly.x PRIVATE POLY_Dmax=${POLY_Dmax})

    palp_sources(CLASS_MAIN_SRC class)
    add_executable(class.x ${CLASS_MAIN_SRC})
    target_link_libraries(class.x PRIVATE objects class_objects)
    target_compile_definitions(class.x PRIVATE POLY_Dmax=${POLY_Dmax})

    palp_sources(CWS_SRC cws LG)
    add_executable(cws.x ${CWS_SRC})
    target_link_libraries(cws.x PRIVATE objects)
    target_compile_definitions(cws.x PRIVATE POLY_Dmax=${POLY_Dmax})

    palp_sources(NEF_MAIN_SRC nef)
    add_executable(nef.x ${NEF_MAIN_SRC})
    target_link_libraries(nef.x PRIVATE objects nef_objects)
    target_compile_definitions(nef.x PRIVATE POLY_Dmax=${POLY_Dmax})

    palp_sources(MORI_MAIN_SRC mori)
    add_executable(mori.x ${MORI_MAIN_SRC})
    target_link_libraries(mori.x PRIVATE objects mori_objects)
    target_compile_definitions(mori.x PRIVATE POLY_Dmax=${POLY_Dmax})

    if (BUILD_TESTING)
        function(add_palp_dimension_targets dim)
            palp_sources(DIM_SOURCES Coord Rat Vertex Polynf)
            add_library(objects-${dim}d STATIC ${DIM_SOURCES})
            target_compile_definitions(objects-${dim}d PUBLIC POLY_Dmax=${dim})
            target_include_directories(objects-${dim}d PUBLIC ${PALP_INCLUDE_DIRS})

            palp_sources(DIM_CLASS_SRC Subpoly Subadd Subdb)
            add_library(class_objects-${dim}d STATIC ${DIM_CLASS_SRC})
            target_compile_definitions(class_objects-${dim}d PUBLIC POLY_Dmax=${dim})
            target_include_directories(class_objects-${dim}d PUBLIC ${PALP_INCLUDE_DIRS})

            palp_sources(DIM_NEF_SRC E_Poly Nefpart LG)
            add_library(nef_objects-${dim}d STATIC ${DIM_NEF_SRC})
            target_compile_definitions(nef_objects-${dim}d PUBLIC POLY_Dmax=${dim})
            target_include_directories(nef_objects-${dim}d PUBLIC ${PALP_INCLUDE_DIRS})

            palp_sources(DIM_MORI_SRC MoriCone SingularInput LG)
            add_library(mori_objects-${dim}d STATIC ${DIM_MORI_SRC})
            target_compile_definitions(mori_objects-${dim}d PUBLIC POLY_Dmax=${dim})
            target_include_directories(mori_objects-${dim}d PUBLIC ${PALP_INCLUDE_DIRS})

            palp_sources(DIM_POLY_SRC poly LG)
            add_executable(poly-${dim}d.x ${DIM_POLY_SRC})
            target_link_libraries(poly-${dim}d.x PRIVATE objects-${dim}d)
            target_compile_definitions(poly-${dim}d.x PRIVATE POLY_Dmax=${dim})

            palp_sources(DIM_CLASS_MAIN_SRC class)
            add_executable(class-${dim}d.x ${DIM_CLASS_MAIN_SRC})
            target_link_libraries(class-${dim}d.x PRIVATE objects-${dim}d class_objects-${dim}d)
            target_compile_definitions(class-${dim}d.x PRIVATE POLY_Dmax=${dim})

            palp_sources(DIM_CWS_SRC cws LG)
            add_executable(cws-${dim}d.x ${DIM_CWS_SRC})
            target_link_libraries(cws-${dim}d.x PRIVATE objects-${dim}d)
            target_compile_definitions(cws-${dim}d.x PRIVATE POLY_Dmax=${dim})

            palp_sources(DIM_NEF_MAIN_SRC nef)
            add_executable(nef-${dim}d.x ${DIM_NEF_MAIN_SRC})
            target_link_libraries(nef-${dim}d.x PRIVATE objects-${dim}d nef_objects-${dim}d)
            target_compile_definitions(nef-${dim}d.x PRIVATE POLY_Dmax=${dim})

            palp_sources(DIM_MORI_MAIN_SRC mori)
            add_executable(mori-${dim}d.x ${DIM_MORI_MAIN_SRC})
            target_link_libraries(mori-${dim}d.x PRIVATE objects-${dim}d mori_objects-${dim}d)
            target_compile_definitions(mori-${dim}d.x PRIVATE POLY_Dmax=${dim})
        endfunction()

        foreach(dim IN LISTS PALP_TEST_DIMENSIONS)
            add_palp_dimension_targets(${dim})
        endforeach()

        add_custom_target(palp_test_tree ALL
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                    ${CMAKE_CURRENT_SOURCE_DIR}/tests
                    ${CMAKE_CURRENT_BINARY_DIR}/tests
            COMMENT "Copying PALP shell tests into the build tree"
        )

        set(PALP_CTEST_ENV)
        if (PALP_TESTS_LONG)
            set(PALP_CTEST_ENV LONG=1)
        endif()

        file(GLOB PALP_TEST_SCRIPTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.sh")
        list(SORT PALP_TEST_SCRIPTS)
        foreach(test_script IN LISTS PALP_TEST_SCRIPTS)
            get_filename_component(test_file "${test_script}" NAME)
            foreach(dim IN LISTS PALP_TEST_DIMENSIONS)
                add_test(
                    NAME "${test_file}-${dim}d"
                    COMMAND ${CMAKE_COMMAND} -E env "DIM=${dim}" ${PALP_CTEST_ENV}
                            "${CMAKE_CURRENT_BINARY_DIR}/tests/${test_file}"
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                )
                set_tests_properties("${test_file}-${dim}d" PROPERTIES DEPENDS palp_test_tree)
            endforeach()
        endforeach()
    endif()
endif()
