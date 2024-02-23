#!/bin/sh
#
# Test the nef -c examples in Section 6.4.11 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 40
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
COMMAND="echo '3 1 1 1 0 0 0 3 0 0 0 1 1 1' | ./nef-${DIM}d.x -f -c3 | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -c3 example on page 40"
EXPECTED=$(cat<<-EOF
3 1 1 1 0 0 0  3 0 0 0 1 1 1 M:100 9 N:7 6  codim=3 #part=7
H:[0] P:0 V0:1 3  V1:4 5
H:[0] P:1 V0:2 3  V1:4 5
np=1 d:1 p:5
EOF
)
if [ $DIM -lt 6 ]; then
   EXPECTED="Please increase POLY_Dmax to at least 6 = 4 + 3 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
fi
run_test


COMMAND="echo '5 1 1 1 1 1' | ./nef-${DIM}d.x -f -c1 | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -c1 example on page 40"
EXPECTED=$(cat<<-EOF
5 1 1 1 1 1 M:126 5 N:6 5  codim=1 #part=1
H:1 101 [-200] P:0
np=1 d:0 p:0
EOF
)
run_test


COMMAND="echo '5 1 1 1 1 1' | ./poly-${DIM}d.x -f"
DESCRIPTION="poly-${DIM}d.x example on page 40"
EXPECTED=$(cat<<-EOF
5 1 1 1 1 1 M:126 5 N:6 5 H:1,101 [-200]
EOF
)
run_test
