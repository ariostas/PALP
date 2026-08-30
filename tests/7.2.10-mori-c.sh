#!/bin/sh
#
# Test the mori -c example in Section 7.2.10 of the PALP manual
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
COMMAND="echo '8 4 1 1 1 1 0 6 3 1 0 1 0 1' | ./mori-${DIM}d.x -fc"
DESCRIPTION="mori-${DIM}d.x -fc example on page 62"
EXPECTED=$(cat<<-EOF
SINGULAR -> divisor classes (integral basis J1 ... J2):
d1=J1+3*J2, d2=J1, d3=-J1+J2, d4=J2, d5=J1, d6=J2
SINGULAR  -> Chern classes of the CY-hypersurface:
c1(CY)=  0
c2(CY)=  10*J1*J2+12*J2^2
c3(CY)=  -252 *[pt]
SINGULAR -> divisor classes (integral basis J1 ... J2):
d1=J1+3*J2, d2=J1, d3=-J1+J2, d4=J2, d5=J1, d6=J2
SINGULAR  -> Chern classes of the CY-hypersurface:
c1(CY)=  0
c2(CY)=  10*J1*J2+12*J2^2
c3(CY)=  -252 *[pt]
EOF
)
run_test "${SKIP}" "${SKIPREASON}"
