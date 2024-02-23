#!/bin/sh
#
# Test the nef -n example in Section 6.4.17 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 46
COMMAND="echo '4 1 1 1 1' | ./nef-${DIM}d.x -f -n"
DESCRIPTION="nef-${DIM}d.x -n example on page 46"
EXPECTED=$(cat<<-EOF
4 1 1 1 1 M:35 4 N:5 4  codim=2 #part=2
3 5  Points of Poly in N-Lattice:
   -1    0    0    1    0
   -1    0    1    0    0
   -1    1    0    0    0
EOF
)
run_test
