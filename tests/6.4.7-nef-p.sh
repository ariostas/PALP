#!/bin/sh
#
# Test the nef -p examples in Section 6.4.7 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 38
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
COMMAND="echo '8 1 1 1 1 1 1 1 1' | ./nef-${DIM}d.x -f -c4 -p | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -c4 -p example on page 38"
EXPECTED=$(cat<<-EOF
8 1 1 1 1 1 1 1 1 M:6435 8 N:9 8  codim=4 #part=5
 P:0 V0:2 3  V1:4 5  V2:6 7
np=1 d:0 p:4
EOF
)
if [ $DIM -lt 10 ]; then
   EXPECTED="Please increase POLY_Dmax to at least 10 = 7 + 4 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
else
    if [ -z "${LONG}" ]; then
	# This test takes forever to run
	SKIP=true
	SKIPREASON="long-running test requires make checklong"
    fi
fi
run_test "${SKIP}" "${SKIPREASON}"
