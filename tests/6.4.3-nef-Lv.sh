#!/bin/sh
#
# Test the nef -Lv examples in Section 6.4.3 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 34
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
COMMAND="./nef-${DIM}d.x -Lv tests/input/6.4.3-nef-Lv.1.txt | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -Lv example on page 34"
EXPECTED=$(cat<<-EOF
M:5 4 N:35 4  codim=2 #part=0
3 4 Vertices in N-lattice:
   -1   -1   -1    3
   -1   -1    3   -1
   -1    3   -1   -1
--------------------
    1    1    1    1  d=4  codim=0
np=0 d:0 p:0
EOF
)
#if [ $DIM -lt 6 ]; then
#   EXPECTED="Please increase POLY_Dmax to at least 6 = 5 + 2 - 1
#(POLY_Dmax >= dim N + codim - 1 is required)"
#fi
run_test



COMMAND="./nef-${DIM}d.x -Lv -N tests/input/6.4.3-nef-Lv.2.txt | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -Lv -N example on pages 34-35"
EXPECTED=$(cat<<-EOF
M:35 4 N:5 4  codim=2 #part=2
3 4 Vertices in N-lattice:
   -1    0    0    1
   -1    0    1    0
   -1    1    0    0
--------------------
    1    1    1    1  d=4  codim=0
H:[0] P:0 V:2 3   (2 2)
np=1 d:0 p:1
EOF
)
run_test
