#!/bin/sh
#
# Test the examples in Section 3.2.11 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 14
COMMAND="echo '5 1 1 1 1 1' | ./poly-${DIM}d.x -fl"
DESCRIPTION="poly-${DIM}d.x -l first example on page 14"
EXPECTED="5 1 1 1 1 1 M:126 5 N:6 5 V:1,101 [-200]"
run_test

COMMAND="echo '3 1 1 1 1 1 1' | ./poly-${DIM}d.x -fl"
DESCRIPTION="poly-${DIM}d.x -l second example on page 14"
EXPECTED="3 1 1 1 1 1 1 M:56 6 F:6 LG: H0:1,0,1 H1:0,20 H2:1 RefI2"
[ $DIM -lt 5 ] && EXPECTED="Increase POLY_Dmax"
run_test

COMMAND="echo '3 1 1 1 1 1 1 /Z3: 0 1 2 0 1 2 3 1 1 1 1 1 1' | ./poly-${DIM}d.x -fl"
DESCRIPTION="poly-${DIM}d.x -l third example on page 14"
EXPECTED="3 1 1 1 1 1 1 /Z3: 0 1 2 0 1 2 M:20 6 F:6 LG: H0:1,0,1 H1:0,20 H2:1 RefI2"
[ $DIM -lt 5 ] && EXPECTED="Increase POLY_Dmax"
run_test
