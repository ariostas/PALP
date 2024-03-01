#!/bin/sh
#
# Test the nef -P examples in Section 6.4.9 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 39
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
COMMAND="echo '4 1 1 1 1' | ./nef-${DIM}d.x -f -P | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -P example on page 39"
EXPECTED=$(cat<<-EOF
4 1 1 1 1 M:35 4 N:5 4  codim=2 #part=2
H:[0] P:0 V:2 3
H:[0] P:1 V:3
np=1 d:0 p:1
EOF
)
run_test
