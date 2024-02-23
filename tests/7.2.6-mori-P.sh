#!/bin/sh
#
# Test the mori -P example in Section 7.2.6 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 60
COMMAND="echo '16 8 4 2 1 1' | ./mori-${DIM}d.x -fP"
DESCRIPTION="mori-${DIM}d.x -fP example on page 60"
EXPECTED=$(cat<<-EOF
4 9  points of P* and IP-simplices
   -1    0    0    2    0    0    0    1    0
   -1    0    0    1    2    0    1    1    0
    0    0    2   -1    1    1    1    0    0
    0    1    1    0   -1    1    0    0    0
-----------------------------------   #IP-simp=3
    8    1    1    4    2    0    0  16=d  codim=0
    4    0    0    2    1    1    0   8=d  codim=1
    2    0    0    1    0    0    1   4=d  codim=2
EOF
)
run_test
