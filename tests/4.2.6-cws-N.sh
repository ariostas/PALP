#!/bin/sh
#
# Test the example in Section 4.2.6 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 24
COMMAND="./cws-${DIM}d.x -N tests/input/4.2.6-cws-N.txt"
DESCRIPTION="cws-${DIM}d.x -N example on page 22"
EXPECTED="5 1 1 1 1 1 /Z5: 4 1 0 0 0 /Z5: 4 0 1 0 0 /Z5: 4 0 0 1 0 "
run_test
