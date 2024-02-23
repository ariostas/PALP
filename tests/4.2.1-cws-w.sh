#!/bin/sh
#
# Test the example in Section 4.2.1 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 22
COMMAND="./cws-${DIM}d.x -w2"
DESCRIPTION="cws-${DIM}d.x -w2 example on page 22"
EXPECTED=$(cat<<-EOF
3  1 1 1  rt
4  1 1 2  rt
6  1 2 3  rt  #=3  #cand=3
EOF
)
run_test
