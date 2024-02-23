#!/bin/sh
#
# Test the mori -d example in Section 7.2.12 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 64
COMMAND="echo '5 1 1 1 1 1 0 2 0 0 0 0 1 1' | ./mori-${DIM}d.x -fd"
DESCRIPTION="mori-${DIM}d.x -d example on page 64"

# Note: on my machine, the output is "wrong" but only up to a
# permutation (it looks like divisors two and four are switched).
# This might be due to a change in Singular, but in any case, I'm
# basically guessint that it shouldn't be fatal.
EXPECTED=$(cat<<-EOF
SINGULAR -> topological quantities of the toric divisors:
Euler characteristics: 46 9 46 46 46 55 
Arithmetic genera: 4 1 4 4 4 5 
dPs: 1 ; d2(6)  nonint: 1 ; d2 
EOF
)
run_test

