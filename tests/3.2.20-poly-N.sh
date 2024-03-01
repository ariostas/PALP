#!/bin/sh
#
# Test the examples in Section 3.2.20 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 16
COMMAND="echo '3402 40 41 486 1134 1701' | ./poly-${DIM}d.x -fN"
DESCRIPTION="poly-${DIM}d.x -N first example on page 16"
# When POLY_Dmax < 6, it looks like it takes fewer permutations
# to arrive at the normal form. This information is purely...
# informational, so it shouldn't trigger a failure.
_perm=43210
[ $DIM -lt 6 ] && _perm=43201
EXPECTED=$(cat<<-EOF
4 5  Normal form of vertices of P    perm=${_perm}
   1   0   0   0 -42
   0   1   0   0 -28
   0   0   1   0 -12
   0   0   0   1  -1
EOF
)
run_test

COMMAND="echo '3486 41 42 498 1162 1743' | ./poly-${DIM}d.x -fN"
DESCRIPTION="poly-${DIM}d.x -N second example on page 16"
# EXPECTED is unchanged from the previous example
run_test
