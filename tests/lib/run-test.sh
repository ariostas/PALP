# The tests compare program output verbatim, and some of them post-process it
# with utilities whose behaviour depends on the locale -- "sort" in particular
# collates punctuation differently outside the C locale, which reorders lines
# that begin with a minus sign. Pin the locale so the expected output in the
# scripts is the only thing that matters.
LC_ALL=C
export LC_ALL

# Succeeds if the named program can be found in PATH. Used by the scripts
# whose examples shell out to an external tool, so that they skip instead of
# failing when the tool is not installed.
have_program() {
    command -v "$1" >/dev/null 2>&1
}

# Run one test and print the result. Most of its arguments
# are passed via environment variables and are required:
#
#   - COMMAND: the test command to run
#   - DESCRIPTION: a short (one-line) description of the test
#   - EXPECTED: the expected output from the test
#
# In addition, this function takes two positional arguments,
#
#   - skip: non-null if the test should be skipped
#   - skipreason: the reason why the test was skipped;
#     required if "skip" is non-null.
#
# These are positional to avoid a common pitfall where e.g. SKIP would
# be set to true at the top of a file, for the first test, and then
# forgotten in subsequent tests causing them to be skipped.
run_test() {
    _SKIP=$1
    _SKIPREASON=$2

    printf "Checking ${DESCRIPTION}... "
    if [ -n "${_SKIP}" ]; then
	printf "SKIP (${_SKIPREASON})"
	printf "\n"
	return
    fi

    # Use eval to support pipes in COMMAND. This isn't
    # a great way to run a stored command, but we also
    # want to _print_ it when the test fails, and that
    # complicates things. If we need to change this in
    # the future, we could instead insist that our caller
    # define (for example) a _test_command() function to
    # perform the test.
    ACTUAL=$(eval "${COMMAND}")

    if [ "${ACTUAL}" = "${EXPECTED}" ]; then
	printf "PASS\n"
    else
	printf "FAIL\n"
	echo "Failed command: ${COMMAND}"
	echo "Expected output"
	echo "---------------"
	echo "${EXPECTED}"
	echo "Actual output"
	echo "-------------"
	echo "${ACTUAL}"
	exit 1
    fi
}
