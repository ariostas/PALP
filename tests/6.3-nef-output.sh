#!/bin/sh
#
# Test the nef example in Section 6.3 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 32
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
INPUT="3 1 1 1 0 0 0 0 0 2 0 0 0 1 1 0 0 0 3 0 0 0 0 0 1 1 1"
COMMAND="echo '${INPUT}' | ./nef-${DIM}d.x -f | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x output example on page 32"
EXPECTED=$(cat<<-EOF
3 1 1 1 0 0 0 0 0  2 0 0 0 1 1 0 0 0  3 0 0 0 0 0 1 1 1 M:300 18 N:9 8  codim=2 #part=15
H:19 19 [0] P:0 V:2 4 6 7
H:9 27 [-36] P:2 V:3 4 6 7
H:3 51 [-96] P:3 V:3 5 6 7
H:3 75 [-144] P:4 V:3 6 7
H:3 51 [-96] P:6 V:4 5 6 7
H:3 51 [-96] P:7 V:4 5 6
H:6 51 [-90] P:8 V:4 6 7
H:3 75 [-144] P:9 V:4 6
H:3 60 [-114] P:10 V:5 6 7
H:3 69 [-132] P:11 V:5 6
H:3 75 [-144] P:12 V:6 7
np=11 d:2 p:2
EOF
)
if [ $DIM -lt 6 ]; then
   EXPECTED="Please increase POLY_Dmax to at least 6 = 5 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
fi
run_test
