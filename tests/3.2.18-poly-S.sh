#!/bin/sh
#
# Test the examples in Section 3.2.18 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 15
COMMAND="echo '5 1 1 1 1 1' | ./poly-${DIM}d.x -fS"
DESCRIPTION="poly-${DIM}d.x -S first example on page 15"
EXPECTED="#GL(Z,4)-Symmetries=120, #VPM-Symmetries=120"
run_test

COMMAND="echo '5 1 1 1 1 1 /Z5: 0 1 2 3 4' | ./poly-${DIM}d.x -fS"
DESCRIPTION="poly-${DIM}d.x -S second example on page 15"
EXPECTED="#GL(Z,4)-Symmetries=20, #VPM-Symmetries=120"
run_test
