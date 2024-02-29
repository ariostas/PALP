#!/bin/sh
#
# Test the nef -m examples in Section 6.4.19 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 48
# The "sed" command strips the unpredictable timing information
COMMAND="echo '14 1 1 1 1 4 6' | ./nef-${DIM}d.x -f -Lv | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -Lv example on page 48"
EXPECTED=$(cat<<-EOF
14 1 1 1 1 4 6 M:1271 13 N:10 8  codim=2 #part=2
5 8 Vertices in N-lattice:
    0   -1    0    0    0    1    0    0
    0   -1    1    0    0    0    0    0
    0   -1    0    1    0    0    0    0
    0   -4    0    0    1    0   -1   -1
    1   -6    0    0    0    0   -1   -2
----------------------------------------
    6    1    1    1    4    1    0    0  d=14  codim=0
    1    0    0    0    1    0    1    0  d=3  codim=3
    2    0    0    0    1    0    0    1  d=4  codim=3
H:1 149 [-296] P:1 V:3 4 5 7  8   (6 8) (1 2) (2 2)
np=1 d:0 p:1
EOF
)
if [ $DIM -lt 6 ]; then
   EXPECTED="Please increase POLY_Dmax to at least 6 = 5 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
fi
run_test

# Pages 48-49
# The "sed" command strips the unpredictable timing information
COMMAND="echo '14 1 1 1 1 4 6 d=2 12' | ./nef-${DIM}d.x -f -Lv -m | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -Lv m example on pages 48-49"
EXPECTED=$(cat<<-EOF
14 1 1 1 1 4 6 d=2 12 M:1270 12 N:11 7  codim=2 #part=2
5 7 Vertices in N-lattice:
    0   -1    0    0    0    1    0
    0   -1    1    0    0    0    0
    0   -1    0    1    0    0    0
    0   -4    0    0    1    0   -2
    1   -6    0    0    0    0   -3
-----------------------------------
    6    1    1    1    4    1    0  d=14  codim=0
    3    0    0    0    2    0    1  d=6  codim=3
 d=12 2H:3 243 [-480] P:0 V:3 5   (2 12) (0 6)
np=1 d:0 p:1
EOF
)
if [ $DIM -lt 6 ]; then
    EXPECTED="Please increase POLY_Dmax (POLY_Dmax >= number of weights is required)"
fi
run_test
