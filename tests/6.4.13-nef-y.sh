#!/bin/sh
#
# Test the nef -y example in Section 6.4.13 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 43
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
COMMAND="./nef-${DIM}d.x -y -N tests/input/6.4.13-nef-y.txt"
DESCRIPTION="nef-${DIM}d.x -y -N example on page 43"
EXPECTED=$(cat<<-EOF
3 4 Vertices of Poly in M-lattice:  M:35 4 N:5 4  codim=2 #part=2
   -1   -1   -1    3
   -1   -1    3   -1
   -1    3   -1   -1
EOF
)
run_test
