cmake_minimum_required(VERSION 3.14...4.1)

project(PALP
  VERSION
    2.21
  DESCRIPTION
    "A Package for Analyzing Lattice Polytopes "
  LANGUAGES
    C
)

set(POLY_Dmax 6 CACHE STRING "Maximum polytope dimension")
set(PALP_TEST_DIMENSIONS "4;5;6;11" CACHE STRING "POLY_Dmax values to build and test with CTest")
set(PALP_TESTS_LONG OFF CACHE BOOL "Run PALP tests that require LONG=1")

# Compiling with -DNDEBUG breaks the code
string(REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCES Coord.c Rat.c Vertex.c Polynf.c)
set(CLASS_SRC Subpoly.c Subadd.c Subdb.c)
set(NEF_SRC E_Poly.c Nefpart.c LG.c)
set(MORI_SRC MoriCone.c SingularInput.c LG.c)

add_library(objects STATIC ${SOURCES})
target_compile_definitions(objects PUBLIC POLY_Dmax=${POLY_Dmax})
target_include_directories(objects PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(class_objects STATIC ${CLASS_SRC})
target_compile_definitions(class_objects PUBLIC POLY_Dmax=${POLY_Dmax})
target_include_directories(class_objects PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(nef_objects STATIC ${NEF_SRC})
target_compile_definitions(nef_objects PUBLIC POLY_Dmax=${POLY_Dmax})
target_include_directories(nef_objects PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(mori_objects STATIC ${MORI_SRC})
target_compile_definitions(mori_objects PUBLIC POLY_Dmax=${POLY_Dmax})
target_include_directories(mori_objects PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if (CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    include(CTest)

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

    if (BUILD_TESTING)
        function(add_palp_dimension_targets dim)
            add_library(objects-${dim}d STATIC ${SOURCES})
            target_compile_definitions(objects-${dim}d PUBLIC POLY_Dmax=${dim})
            target_include_directories(objects-${dim}d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

            add_library(class_objects-${dim}d STATIC ${CLASS_SRC})
            target_compile_definitions(class_objects-${dim}d PUBLIC POLY_Dmax=${dim})
            target_include_directories(class_objects-${dim}d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

            add_library(nef_objects-${dim}d STATIC ${NEF_SRC})
            target_compile_definitions(nef_objects-${dim}d PUBLIC POLY_Dmax=${dim})
            target_include_directories(nef_objects-${dim}d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

            add_library(mori_objects-${dim}d STATIC ${MORI_SRC})
            target_compile_definitions(mori_objects-${dim}d PUBLIC POLY_Dmax=${dim})
            target_include_directories(mori_objects-${dim}d PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

            add_executable(poly-${dim}d.x poly.c LG.c)
            target_link_libraries(poly-${dim}d.x PRIVATE objects-${dim}d)
            target_compile_definitions(poly-${dim}d.x PRIVATE POLY_Dmax=${dim})

            add_executable(class-${dim}d.x class.c)
            target_link_libraries(class-${dim}d.x PRIVATE objects-${dim}d class_objects-${dim}d)
            target_compile_definitions(class-${dim}d.x PRIVATE POLY_Dmax=${dim})

            add_executable(cws-${dim}d.x cws.c LG.c)
            target_link_libraries(cws-${dim}d.x PRIVATE objects-${dim}d)
            target_compile_definitions(cws-${dim}d.x PRIVATE POLY_Dmax=${dim})

            add_executable(nef-${dim}d.x nef.c)
            target_link_libraries(nef-${dim}d.x PRIVATE objects-${dim}d nef_objects-${dim}d)
            target_compile_definitions(nef-${dim}d.x PRIVATE POLY_Dmax=${dim})

            add_executable(mori-${dim}d.x mori.c)
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
            endforeach()
        endforeach()
    endif()
endif()
