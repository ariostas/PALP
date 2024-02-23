#!/bin/sh
#
# Test the nef -V example in Section 6.4.21 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 49-50
# The "sed" command strips the unpredictable timing information
COMMAND="echo '4 1 1 1 1' | ./nef-${DIM}d.x -f -V | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -V example on pages 49-50"
EXPECTED=$(cat<<-EOF
4 1 1 1 1 M:35 4 N:5 4  codim=2 #part=2
3 4  Vertices of P:
   -1    0    0    1
   -1    0    1    0
   -1    1    0    0
H:[0] P:0 V:2 3
np=1 d:0 p:1
EOF
)
run_test
