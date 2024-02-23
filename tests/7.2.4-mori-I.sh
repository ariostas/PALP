#!/bin/sh
#
# Test the mori -I example in Section 7.2.4 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 59
COMMAND="echo '8 4 1 1 1 1 0 6 3 1 0 1 0 1' | ./mori-${DIM}d.x -fI"
DESCRIPTION="mori-${DIM}d.x -fI example on page 59"
EXPECTED=$(cat<<-EOF
Incidence: 110101 111100 011111 101011 101110 100111 111001
EOF
)
run_test
