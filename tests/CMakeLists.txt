# CTest wiring for the shell test suite, mirroring the "check" and "checklong"
# targets of the GNUmakefile.
#
# The scripts in this directory follow the protocol documented in tests/README:
# they source tests/lib/run-test.sh, read the input files under tests/input,
# and invoke ./<program>-<DIM>d.x -- all relative to the current directory.
# We therefore build the dimension-specific programs into a dedicated
# directory, link the tests/ tree next to them, and run each script from
# there. The scripts exit non-zero on the first failing check, which is
# exactly what CTest wants.

set(PALP_TEST_DIMENSIONS "4;5;6;11" CACHE STRING
  "POLY_Dmax values to build test programs for")

if (NOT PALP_TEST_DIMENSIONS)
  message(FATAL_ERROR "PALP_TEST_DIMENSIONS must not be empty")
endif()
foreach(dim IN LISTS PALP_TEST_DIMENSIONS)
  if (NOT "${dim}" MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR
      "PALP_TEST_DIMENSIONS must contain positive integers, got '${dim}'")
  endif()
endforeach()

find_program(PALP_SH_EXECUTABLE NAMES sh)
if (NOT PALP_SH_EXECUTABLE)
  message(FATAL_ERROR
    "PALP_BUILD_TESTS requires a POSIX shell, but no 'sh' was found")
endif()

# Directory the test scripts run in: it holds the *-<dim>d.x programs and a
# link to this source directory, so that both "./poly-6d.x" and
# "tests/input/..." resolve.
set(PALP_TEST_WORKDIR "${PROJECT_BINARY_DIR}/test-tree")
file(MAKE_DIRECTORY "${PALP_TEST_WORKDIR}")
file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}" "${PALP_TEST_WORKDIR}/tests"
     SYMBOLIC RESULT palp_link_result)
if (NOT palp_link_result STREQUAL "0")
  # Platforms without usable symlinks: fall back to a copy. The inputs are
  # static, so this only needs redoing when CMake re-runs.
  file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/" DESTINATION "${PALP_TEST_WORKDIR}/tests")
endif()

foreach(dim IN LISTS PALP_TEST_DIMENSIONS)
  palp_add_libraries("-${dim}d" ${dim})
  palp_add_programs("-${dim}d")

  # Force the programs into PALP_TEST_WORKDIR for every generator, including
  # multi-config ones which would otherwise add a per-configuration subdir.
  set_target_properties(${PALP_PROGRAM_TARGETS} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${PALP_TEST_WORKDIR}")
  foreach(config IN LISTS CMAKE_CONFIGURATION_TYPES)
    string(TOUPPER "${config}" config_upper)
    set_target_properties(${PALP_PROGRAM_TARGETS} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY_${config_upper} "${PALP_TEST_WORKDIR}")
  endforeach()

  list(APPEND PALP_ALL_TEST_PROGRAMS ${PALP_PROGRAM_TARGETS})
endforeach()

file(GLOB PALP_TEST_SCRIPTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.sh")
list(SORT PALP_TEST_SCRIPTS)

foreach(script IN LISTS PALP_TEST_SCRIPTS)
  # NAME_WE would cut at the first dot, and every script name contains
  # several (e.g. "6.4.4-nef-H.sh"), so strip just the extension.
  get_filename_component(name "${script}" NAME)
  string(REGEX REPLACE "\\.sh$" "" name "${name}")

  # A handful of scripts gate their slowest checks on LONG being set; those
  # get a second, separately labelled test so that "checklong" can pick them
  # up while a plain "check" stays fast.
  file(READ "${script}" script_contents)
  set(script_modes "")
  if (script_contents MATCHES "\\\$\\{?LONG")
    set(script_modes long)
  endif()

  foreach(dim IN LISTS PALP_TEST_DIMENSIONS)
    add_test(NAME "${name}-${dim}d"
      COMMAND "${CMAKE_COMMAND}" -E env "DIM=${dim}"
              "${PALP_SH_EXECUTABLE}" "${script}"
      WORKING_DIRECTORY "${PALP_TEST_WORKDIR}")
    set_tests_properties("${name}-${dim}d" PROPERTIES
      LABELS "palp;dim${dim}")

    foreach(mode IN LISTS script_modes)
      add_test(NAME "${name}-${dim}d-${mode}"
        COMMAND "${CMAKE_COMMAND}" -E env "DIM=${dim}" "LONG=1"
                "${PALP_SH_EXECUTABLE}" "${script}"
        WORKING_DIRECTORY "${PALP_TEST_WORKDIR}")
      # No default 1500s timeout: these are the checks the suite deliberately
      # skips unless asked for, and they can run for hours.
      set_tests_properties("${name}-${dim}d-${mode}" PROPERTIES
        LABELS "palp;dim${dim};long" TIMEOUT 0)
    endforeach()
  endforeach()
endforeach()

list(LENGTH PALP_TEST_SCRIPTS PALP_TEST_SCRIPT_COUNT)
message(STATUS
  "PALP tests: ${PALP_TEST_SCRIPT_COUNT} scripts x dimensions [${PALP_TEST_DIMENSIONS}]")

# Convenience targets matching "make check" / "make checklong".
add_custom_target(check
  COMMAND "${CMAKE_CTEST_COMMAND}" --build-config $<CONFIG>
          --output-on-failure --label-exclude long
  WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
  COMMENT "Running the PALP test suite"
  USES_TERMINAL)
add_dependencies(check ${PALP_ALL_TEST_PROGRAMS})

add_custom_target(checklong
  COMMAND "${CMAKE_CTEST_COMMAND}" --build-config $<CONFIG>
          --output-on-failure
  WORKING_DIRECTORY "${PROJECT_BINARY_DIR}"
  COMMENT "Running the PALP test suite, including the long-running checks"
  USES_TERMINAL)
add_dependencies(checklong ${PALP_ALL_TEST_PROGRAMS})
