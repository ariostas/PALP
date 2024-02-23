#!/bin/sh
#
# Test the examples in Section 2.2 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 6, only makes sense for POLY_Dmax < 7.
if [ $DIM -ge 7 ]; then
     SKIP=true
     SKIPREASON="POLY_Dmax too large"
fi
COMMAND="echo '8 1 1 1 1 1 1 1 1' | ./poly-${DIM}d.x -f"
DESCRIPTION="poly-${DIM}d.x POLY_Dmax example on page 6"
EXPECTED="Please increase POLY_Dmax to at least 7"
run_test "${SKIP}" "${SKIPREASON}"

# Page 8, nef will fail and complain about POLY_Dmax for all
# current values of DIM. We can't "fix" this example because
# we would need to change VERT_Nmax as well.
COMMAND="echo '7 9' | ./nef-${DIM}d.x -f -Lp -N -c6 -P"
DESCRIPTION="nef-${DIM}d.x example on page 8"
EXPECTED=$(cat<<-EOF
Please increase POLY_Dmax to at least 12 = 7 + 6 - 1
(POLY_Dmax >= dim N + codim - 1 is required)
EOF
)
run_test
