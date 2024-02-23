#!/bin/sh
#
# Test the nef -D examples in Section 6.4.8 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 38
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
INPUT="3 1 1 1 0 0 0 3 0 0 0 1 1 1"
COMMAND="echo '${INPUT}' | ./nef-${DIM}d.x -f -D | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -D example on page 38"
EXPECTED=$(cat<<-EOF
3 1 1 1 0 0 0  3 0 0 0 1 1 1 M:100 9 N:7 6  codim=2 #part=5
H:4 [0] h1=2 P:0 V:2 3 5   D
H:20 [24] P:1 V:3 4 5
H:20 [24] P:2 V:3 5
H:20 [24] P:3 V:4 5
np=3 d:1 p:1
EOF
)
if [ $DIM -lt 5 ]; then
   EXPECTED="Please increase POLY_Dmax to at least 5 = 4 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
fi
run_test
