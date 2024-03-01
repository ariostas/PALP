#!/bin/sh
#
# Test the nef -R example in Section 6.4.20 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 49
INPUT="6 3 2 1 0 0 5 0 0 1 1 3"
COMMAND="echo '${INPUT}' | ./nef-${DIM}d.x -f -R"
DESCRIPTION="nef-${DIM}d.x -R example on page 49"
EXPECTED=$(cat<<-EOF
3 7  Vertices of input polytope:
   -1    1    0    0    1    0   -1
    0    1   -1    1    4    4    1
   -1    0    0    0   -1   -1   -1
EOF
)
run_test
