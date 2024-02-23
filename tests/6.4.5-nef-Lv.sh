#!/bin/sh
#
# Test the nef -Lv examples in Section 6.4.5 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 37
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
INPUT="5 1 1 1 1 1 0 0 4 0 0 0 1 1 1 1"
COMMAND="echo '${INPUT}' | ./nef-${DIM}d.x -f -Lv | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -Lv example on page 37"
EXPECTED=$(cat<<-EOF
5 1 1 1 1 1 0 0  4 0 0 0 1 1 1 1 M:378 12 N:8 7  codim=2 #part=8
5 7 Vertices in N-lattice:
    0   -1    0    1    0    0    0
    0   -1    1    0    0    0    0
   -1    0    0    0    0    0    1
   -1    1    0    0    1    0    0
   -1    1    0    0    0    1    0
-----------------------------------
    1    1    1    1    0    0    1  d=5  codim=1
    1    0    0    0    1    1    1  d=4  codim=2
H:2 64 [-124] P:0 V:0 6   (2 3) (2 2)
H:2 64 [-124] P:1 V:0 1 6   (3 2) (2 2)
H:2 74 [-144] P:2 V:2 3 5   (2 3) (1 3)
H:2 64 [-124] P:3 V:3 5 6   (2 3) (2 2)
H:2 86 [-168] P:4 V:3 5   (1 4) (1 3)
H:2 74 [-144] P:5 V:3 6   (2 3) (1 3)
np=6 d:0 p:2
EOF
)
if [ $DIM -lt 6 ]; then
   EXPECTED="Please increase POLY_Dmax to at least 6 = 5 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
fi
run_test
