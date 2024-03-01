#!/bin/sh
#
# Test the example in Section 3.2.23 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 17
COMMAND="echo '6 1 2 3' | ./poly-${DIM}d.x -fP"
DESCRIPTION="poly-${DIM}d.x -P example on page 17"
EXPECTED=$(cat<<-EOF
2 7  points of P-dual and IP-simplices
    1    0   -2   -1    0   -1    0
    0    1   -3   -2   -1   -1    0
------------------------------    #IP-simp=4
    2    3    1    0    0    0   6=d  codim=0
    1    2    0    1    0    0   4=d  codim=0
    1    1    0    0    0    1   3=d  codim=0
    0    1    0    0    1    0   2=d  codim=1
EOF
)
if [ $DIM -lt 6 ]; then
    # In this example, the first and third COLUMNS
    # are switched for POLY_Dmax < 6.
    EXPECTED=$(cat<<-EOF
2 7  points of P-dual and IP-simplices
   -2    0    1   -1    0   -1    0
   -3    1    0   -2   -1   -1    0
------------------------------    #IP-simp=4
    1    3    2    0    0    0   6=d  codim=0
    0    2    1    1    0    0   4=d  codim=0
    0    1    1    0    0    1   3=d  codim=0
    0    1    0    0    1    0   2=d  codim=1
EOF
   )
fi
run_test
