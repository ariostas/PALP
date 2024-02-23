#!/bin/sh
#
# Test the nef -Lp examples in Section 6.4.6 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 37-38
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
INPUT="5 1 1 1 1 1 0 0 10 2 2 2 2 0 1 1"
COMMAND="echo '${INPUT}' | ./nef-${DIM}d.x -f -Lp | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -Lp example on pages 37-38"
EXPECTED=$(cat<<-EOF
5 1 1 1 1 1 0 0  10 2 2 2 2 0 1 1 M:378 6 N:8 6  codim=2 #part=4
5 8  Points of Poly in N-Lattice:
   -1    0    0    0    1    0    0    0
   -1    0    1    0    0    0    0    0
   -1    0    0    1    0    0    0    0
   -1    2    0    0    0    0    1    0
   -1    1    0    0    0    1    1    0
----------------------------------------
    2    1    2    2    2    1    0  d=10  codim=0
    1    0    1    1    1    0    1  d=5  codim=1
H:2 86 [-168] P:0 V:1 5  6   (2 8) (1 4)
H:2 68 [-132] P:1 V:2 3 4   (6 4) (3 2)
H:2 68 [-132] P:2 V:3 4   (4 6) (2 3)
np=3 d:0 p:1
EOF
)
if [ $DIM -lt 6 ]; then
   EXPECTED="Please increase POLY_Dmax to at least 6 = 5 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
fi
run_test
