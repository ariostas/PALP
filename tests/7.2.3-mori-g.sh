#!/bin/sh
#
# Test the mori -g example in Section 7.2.3 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 58
COMMAND="echo '8 4 1 1 1 1 0 6 3 1 0 1 0 1' | ./mori-${DIM}d.x -fg"
DESCRIPTION="mori-${DIM}d.x -fg example on page 58"
EXPECTED=$(cat<<-EOF
8 Triangulation
110101 111100 101011 101110 100111 111001 001111 011101
2 SR-ideal
010010 101101
9 Triangulation
110101 111100 101011 101110 100111 111001 010111 011011 011110
2 SR-ideal
110010 001101
EOF
)
run_test
