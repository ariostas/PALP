#!/bin/sh
#
# Test the example in Section 3.2.27 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 19-20
# We sort the output because for lower POLY_Dmax, it's
# violently permuted. The "-r" simply ensures that the
# vol=... header appears before the coordinates.
COMMAND="./poly-${DIM}d.x -B2 tests/input/3.2.27-poly-B.txt | sort -br"
DESCRIPTION="poly-${DIM}d.x -B2 example on pages 19-20"
EXPECTED=$(cat<<-EOF
vol=5, baricent=(5,0,0,0,0)/6
IPs:
 2 2 0 0 0  cd=4
 2 1 1 0 0  cd=3
 2 1 0 1 0  cd=3
 2 1 0 0 1  cd=3
 2 1 0 0 0  cd=0
 2 0 2 0 0  cd=4
 2 0 1 1 0  cd=3
 2 0 1 0 1  cd=3
 2 0 1 0 0  cd=0
 2 0 0 2 0  cd=4
 2 0 0 1 1  cd=3
 2 0 0 1 0  cd=0
 2 0 0 0 2  cd=4
 2 0 0 0 1  cd=0
 2 0 0 0 0  cd=0
 2 0 -1 -1 -1  cd=3
 2 -2 -2 -2 -2  cd=4
 2 -1 0 -1 -1  cd=3
 2 -1 -1 0 -1  cd=3
 2 -1 -1 -1 0  cd=3
 2 -1 -1 -1 -1  cd=0
 1 1 0 0 0  cd=4
 1 0 1 0 0  cd=4
 1 0 0 1 0  cd=4
 1 0 0 0 1  cd=4
 1 0 0 0 0  cd=0
 1 -1 -1 -1 -1  cd=4
 0 0 0 0 0  cd=5
EOF
)
[ $DIM -lt 5 ] && EXPECTED="Please increase POLY_Dmax to at least 5"
run_test
