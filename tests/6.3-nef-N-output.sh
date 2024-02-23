#!/bin/sh
#
# Test the nef -N example in Section 6.3 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 33
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
COMMAND="./nef-${DIM}d.x -N tests/input/6.3-nef-N-output.txt | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -N output example on page 33"
EXPECTED=$(cat<<-EOF
M:300 18 N:9 8  codim=2 #part=15
H:3 51 [-96] P:0 V:2 3 4 7
H:3 51 [-96] P:1 V:2 4 6 7
H:3 60 [-114] P:2 V:2 4 7
H:3 51 [-96] P:3 V:2 6 7
H:3 69 [-132] P:4 V:2 7
H:9 27 [-36] P:5 V:3 4 6 7
H:3 75 [-144] P:6 V:3 4 7
H:19 19 [0] P:8 V:4 5 6 7
H:6 51 [-90] P:9 V:4 6 7
H:3 75 [-144] P:10 V:4 7
H:3 75 [-144] P:13 V:6 7
np=11 d:2 p:2
EOF
)
if [ $DIM -lt 6 ]; then
   EXPECTED="Please increase POLY_Dmax to at least 6 = 5 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
fi
run_test
