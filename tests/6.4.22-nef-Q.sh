#!/bin/sh
#
# Test the nef -Q examples in Section 6.4.22 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 50
# The "sed" command strips the unpredictable timing information
INPUT="3 1 1 1 0 0 0 3 0 0 0 1 1 1"
COMMAND="echo '${INPUT}' | ./nef-${DIM}d.x -f -Q | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -Q first example on page 50"
EXPECTED=$(cat<<-EOF
3 1 1 1 0 0 0  3 0 0 0 1 1 1 M:100 9 N:7 6  codim=2 #part=5
H:4 [0] h1=2 P:0 V:2 3 5   D
np=4 d:1 p:0
EOF
)
if [ $DIM -lt 5 ]; then
    EXPECTED="Please increase POLY_Dmax to at least 5 = 4 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
fi
run_test

# Page 50
# The "sed" command strips the unpredictable timing information
COMMAND="echo '4 1 1 1 1' | ./nef-${DIM}d.x -f -Q | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -Q second example on page 50"
EXPECTED=$(cat<<-EOF
4 1 1 1 1 M:35 4 N:5 4  codim=2 #part=2
np=2 d:0 p:0
EOF
)
run_test
