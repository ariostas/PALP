#!/bin/sh
#
# Test the mori -b example in Section 7.2.8 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Every example below makes mori shell out to Singular ("Singular -q < ...",
# see SingularInput.c), which is an optional external dependency.
if ! have_program Singular; then
    SKIP=true
    SKIPREASON="Singular is not installed"
fi

# Page 61
COMMAND="echo '4 1 1 1 1' | ./mori-${DIM}d.x -bf"
DESCRIPTION="mori-${DIM}d.x -bf first example on page 61"
EXPECTED=$(cat<<-EOF
SINGULAR  -> Arithmetic genera and Euler number of the CY:
chi_0:  2 , chi_1: -20  [ 24 ]
EOF
)
run_test "${SKIP}" "${SKIPREASON}"

COMMAND="echo '8 4 1 1 1 1 0 6 3 1 0 1 0 1' | ./mori-${DIM}d.x -bf"
DESCRIPTION="mori-${DIM}d.x -bf second example on page 61"
EXPECTED=$(cat<<-EOF
SINGULAR  -> Arithmetic genera and Euler number of the CY:
chi_0:  0 , chi_1: 126  [ -252 ]
SINGULAR  -> Arithmetic genera and Euler number of the CY:
chi_0:  0 , chi_1: 126  [ -252 ]
EOF
)
run_test "${SKIP}" "${SKIPREASON}"
