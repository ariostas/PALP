#!/bin/sh
#
# Test the mori -D example in Section 7.2.14 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 65
COMMAND="./mori-${DIM}d.x -DP tests/input/7.2.14-mori-D.txt"
DESCRIPTION="mori-${DIM}d.x -DP example on page 65"
EXPECTED=$(cat<<-EOF
4 9  points of P* and IP-simplices
   -1    2    0    0    0    0    0    1    0
   -1    1    2    0    0    0    1    1    0
    0   -1    1    0    2    1    1    0    0
    0    0   -1    1    1    1    0    0    0
-----------------------------------   #IP-simp=3
    8    4    2    1    1    0    0  16=d  codim=0
    4    2    1    0    0    1    0   8=d  codim=1
    2    1    0    0    0    0    1   4=d  codim=2
EOF
)
run_test

