#!/bin/sh
#
# Test the example in Section 3.2.25 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 18. The expected output in the manual is truncated, but the
# part that's there should agree.
COMMAND="echo '84 1 1 12 28 42' | ./poly-${DIM}d.x -f -12"
DESCRIPTION="poly-${DIM}d.x -12 example on pages 18"
EXPECTED=$(cat<<-EOF
4 25  Em:7 3 n:7 3  Km:24 4 n:24 4  M:680 5 N:26 5  p=13bgjn256789acdefhiklmo04
 1  0 -2 -1  0 -1  0 -14 -12 -10 -8 -6 -4 -9 -7 -5 -3 -4 -2 -7 -5 -3 -2  0 -28
 0  1 -3 -2 -1 -1  0 -21 -18 -15 -12 -9 -6 -14 -11 -8 -5 -7 -4 -10 -7 -4 -3  0 -42
 0  0  0  0  0  0  1 -6 -5 -4 -3 -2 -1 -4 -3 -2 -1 -2 -1 -3 -2 -1 -1  0 -12
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1
EOF
)
if [ $DIM -lt 6 ]; then
    # When POLY_Dmax < 6, the first two rows are swapped, but only
    # after the first two columns. This is a little weird to me but
    # it's not impossible if e.g. the first two represent a relabeling
    # of a basis in a subspace orthogonal to the rest.
    EXPECTED=$(cat<<-EOF
4 25  Em:7 3 n:7 3  Km:24 4 n:24 4  M:680 5 N:26 5  p=13bgjn256789acdefhiklmo04
 1  0 -3 -2 -1 -1  0 -21 -18 -15 -12 -9 -6 -14 -11 -8 -5 -7 -4 -10 -7 -4 -3  0 -42
 0  1 -2 -1  0 -1  0 -14 -12 -10 -8 -6 -4 -9 -7 -5 -3 -4 -2 -7 -5 -3 -2  0 -28
 0  0  0  0  0  0  1 -6 -5 -4 -3 -2 -1 -4 -3 -2 -1 -2 -1 -3 -2 -1 -1  0 -12
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1
EOF
    )
fi
run_test
