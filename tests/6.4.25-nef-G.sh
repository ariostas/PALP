#!/bin/sh
#
# Test the nef -G example in Section 6.4.25 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}


# Page 54
COMMAND="./nef-${DIM}d.x -G tests/input/6.4.25-nef-G.1.txt"
DESCRIPTION="nef-${DIM}d.x -G first example on page 54"
EXPECTED=$(cat<<-EOF
M:4 4 N:4 4 H:[0] h0=0
EOF
)
run_test

COMMAND="echo '3 1 1 1 1 1 1' | ./nef-${DIM}d.x -f -G"
DESCRIPTION="nef-${DIM}d.x -G second example on page 54"
EXPECTED=$(cat<<-EOF
3 1 1 1 1 1 1 M:56 6 N:6 6 H:20 [24]
EOF
)
if [ $DIM -lt 6 ]; then
    EXPECTED="Please increase POLY_Dmax to at least 6 = 5 + 1
(option -G requires POLY_Dmax >= dim(cone) = dim(support) + 1)"
fi
run_test


COMMAND="echo '7 1 1 1 2 3 3 3' | ./nef-${DIM}d.x -f -G"
DESCRIPTION="nef-${DIM}d.x -G third example on page 54"
EXPECTED=$(cat<<-EOF
7 1 1 1 2 3 3 3 M:154 18 F:9 
EOF
)
if [ $DIM -lt 7 ]; then
    EXPECTED="Please increase POLY_Dmax to at least 7 = 6 + 1
(option -G requires POLY_Dmax >= dim(cone) = dim(support) + 1)"
fi
run_test


COMMAND="echo '7 1 1 2 2 2 3 3' | ./nef-${DIM}d.x -f -G"
DESCRIPTION="nef-${DIM}d.x -G fourth example on page 54"
EXPECTED=$(cat<<-EOF
7 1 1 2 2 2 3 3 M:116 18 N:9 9 H:2 70 [-136]
EOF
)
if [ $DIM -lt 7 ]; then
    EXPECTED="Please increase POLY_Dmax to at least 7 = 6 + 1
(option -G requires POLY_Dmax >= dim(cone) = dim(support) + 1)"
fi
run_test


COMMAND="./nef-${DIM}d.x -G tests/input/6.4.25-nef-G.2.txt"
DESCRIPTION="nef-${DIM}d.x -G fifth example on page 54"
EXPECTED=$(cat<<-EOF
Warning: Input has index 3, should be 2!   M:3 3 F:3 
EOF
)
run_test

INPUT="1 1 1 0 0 0 0 2 0 0 1 1 1 1"
COMMAND="echo '${INPUT}' | ./nef-${DIM}d.x -f -G"
DESCRIPTION="nef-${DIM}d.x -G sixth example on page 55"
EXPECTED=$(cat<<-EOF
1 1 1 0 0 0 0  2 0 0 1 1 1 1 M:20 8 N:6 6 H:[0]
EOF
)
if [ $DIM -lt 5 ]; then
    EXPECTED="Please increase POLY_Dmax to at least 5 = 4 + 1
(option -G requires POLY_Dmax >= dim(cone) = dim(support) + 1)"
fi
run_test
