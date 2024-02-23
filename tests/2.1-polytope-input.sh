#!/bin/sh
#
# Test the examples in Section 2.1 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 5: Test both the transposed and non-transposed inputs.
for t in "" "-transpose"; do
    COMMAND="./poly-${DIM}d.x tests/input/2.1-polytope-input${t}.txt"
    DESCRIPTION="poly-${DIM}d.x example on page 5"
    [ -n "${t}" ] && DESCRIPTION="${DESCRIPTION} (transposed)"
    EXPECTED="M:6 3 F:3"
    run_test
done

# Page 5/6
COMMAND="echo '5 1 1 1 1 1' | ./poly-${DIM}d.x -fv tests/input/2.1-polytope-input-b.txt"
DESCRIPTION="poly-${DIM}d.x -v example on pages 5/6"
EXPECTED=$(cat<<-EOF
4 5  Vertices of P
   -1    4   -1   -1   -1
   -1   -1    4   -1   -1
   -1   -1   -1    4   -1
   -1   -1   -1   -1    4
EOF
)
run_test

# Page 6
COMMAND="echo '2 1 1 0 0 2 0 0 1 1' | ./poly-${DIM}d.x -fv tests/input/2.1-polytope-input-c.txt"
DESCRIPTION="poly-${DIM}d.x -v first example on page 6"
EXPECTED=$(cat<<-EOF
2 4  Vertices of P
   -1    1   -1    1
   -1   -1    1    1
EOF
)
run_test


COMMAND="echo '2 1 1 0 0 2 0 0 1 1 /Z2: 1 0 1 0' | ./poly-${DIM}d.x -fv tests/input/2.1-polytope-input-d.txt"
DESCRIPTION="poly-${DIM}d.x -v second example on page 6"
EXPECTED=$(cat<<-EOF
2 4  Vertices of P
   -1    0    0    1
    1   -1    1   -1
EOF
)
run_test
