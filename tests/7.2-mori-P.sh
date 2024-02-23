#!/bin/sh
#
# Test the mori -P example in Section 7.2 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 57-58
INPUT="8 4 1 1 1 1 0 6 3 1 0 1 0 1"
COMMAND="echo '${INPUT}' | ./mori-${DIM}d.x -fP"
DESCRIPTION="mori-${DIM}d.x -P example on pages 57-58"
EXPECTED=$(cat<<-EOF
4 8  points of P* and IP-simplices
   -1    0    0    0    1    3    1    0
    0    0    0    1    0   -1    0    0
   -1    1    0    0    0    3    1    0
    1    0    1    0    0   -4   -1    0
------------------------------   #IP-simp=2
    4    1    0    1    1    1   8=d  codim=0
    3    0    1    1    0    1   6=d  codim=1
EOF
)
run_test

# Might as well include the repeated example from section 7.2.2, which
# is the same as above except using -f and a pipeline rather than an
# input file.
COMMAND="echo '8 4 1 1 1 1 0 6 3 1 0 1 0 1' | ./mori-${DIM}d.x -fP"
DESCRIPTION="mori-${DIM}d.x -fP example on page 58"
run_test
