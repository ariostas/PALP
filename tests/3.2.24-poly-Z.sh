#!/bin/sh
#
# Test the example in Section 3.2.24 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 17-18
COMMAND="./poly-${DIM}d.x -VZD tests/input/3.2.24-poly-Z.txt"
DESCRIPTION="poly-${DIM}d.x -Z example on pages 17-18"
EXPECTED=$(cat<<-EOF
3 5  vertices of P-dual and IP-simplices
   -1   -1    2    0    0
   -1    2   -1    0    0
    0    0    0    1   -1
-------------------------   #IP-simp=2 I=3 /Z3: 2 1 0 0 0
    1    1    1    0    0    3=d  codim=1 /Z3: 2 1 0 0 0
    0    0    0    1    1    2=d  codim=2
EOF
)
run_test
