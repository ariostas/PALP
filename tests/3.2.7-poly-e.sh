#!/bin/sh
#
# Test the examples in Section 3.2.7 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 13
COMMAND="./poly-${DIM}d.x -e tests/input/3.2.7-poly-e.1.txt"
DESCRIPTION="poly-${DIM}d.x -e first example on page 13"
EXPECTED=$(cat<<-EOF
3 2  Equations of P
   1   0     0
   0   1     0
  -1  -1     1
EOF
)
if [ $DIM -lt 6 ]; then
    # The output equations are permuted when POLY_Dmax < 6 (the first
    # and third rows are swapped), but they're the same equations.
    EXPECTED=$(cat<<-EOF
3 2  Equations of P
  -1  -1     1
   0   1     0
   1   0     0
EOF
    )
fi
run_test

COMMAND="./poly-${DIM}d.x -e tests/input/3.2.7-poly-e.2.txt"
DESCRIPTION="poly-${DIM}d.x -e second example on page 13"
EXPECTED=$(cat<<-EOF
3 2  Vertices of P-dual <-> Equations of P
   2  -1
  -1   2
  -1  -1
EOF
)
if [ $DIM -lt 6 ]; then
    # The output equations are permuted when POLY_Dmax < 6 (the first
    # and third rows are swapped), but they're the same equations.
    EXPECTED=$(cat<<-EOF
3 2  Vertices of P-dual <-> Equations of P
  -1  -1
  -1   2
   2  -1
EOF
    )
fi
run_test
