#!/bin/sh
#
# Test the mori -i example in Section 7.2.9 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Every example below makes mori shell out to Singular ("Singular -q < ...",
# see SingularInput.c), which is an optional external dependency.
if ! have_program Singular; then
    SKIP=true
    SKIPREASON="Singular is not installed"
fi

# Page 62
COMMAND="echo '8 4 1 1 1 1 0 6 3 1 0 1 0 1' | ./mori-${DIM}d.x -fi"
DESCRIPTION="mori-${DIM}d.x -fi example on page 62"
EXPECTED=$(cat<<-EOF
SINGULAR -> divisor classes (integral basis J1 ... J2):
d1=J1+3*J2, d2=J1, d3=-J1+J2, d4=J2, d5=J1, d6=J2
SINGULAR -> intersection polynomial:
2*J1*J2^2+2*J2^3
SINGULAR -> divisor classes (integral basis J1 ... J2):
d1=J1+3*J2, d2=J1, d3=-J1+J2, d4=J2, d5=J1, d6=J2
SINGULAR -> intersection polynomial:
2*J1*J2^2+2*J2^3
EOF
)
run_test "${SKIP}" "${SKIPREASON}"
