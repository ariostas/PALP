#!/bin/sh
#
# Test the nef -v examples in Section 6.4.18 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 46
COMMAND="./nef-${DIM}d.x -v tests/input/6.4.18-nef-v.txt"
DESCRIPTION="nef-${DIM}d.x -v example on page 46"
EXPECTED=$(cat<<-EOF
3 4 P:35 E   -1    3   -1   -1E   -1   -1    3   -1E   -1   -1   -1    3
4 9 P:100 E   -1    2   -1   -1    2   -1   -1    2   -1E   -1   -1    2   -1   -1    2   -1   -1    2E   -1   -1   -1    2    2    2   -1   -1   -1E   -1   -1   -1   -1   -1   -1    2    2    2


2  of  2

  35#    1
 100#    1
EOF
)
if [ $DIM -lt 5 ]; then
   EXPECTED="3 4 P:35 E   -1    3   -1   -1E   -1   -1    3   -1E   -1   -1   -1    3
Please increase POLY_Dmax to at least 5 = 4 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
fi
run_test

# Page 46
COMMAND="./nef-${DIM}d.x -v -u50 tests/input/6.4.18-nef-v.txt"
DESCRIPTION="nef-${DIM}d.x -v -u50 example on page 46"
EXPECTED=$(cat<<-EOF
3 4 P:35 E   -1    3   -1   -1E   -1   -1    3   -1E   -1   -1   -1    3


1  of  2

  35#    1
EOF
)
if [ $DIM -lt 5 ]; then
   EXPECTED="3 4 P:35 E   -1    3   -1   -1E   -1   -1    3   -1E   -1   -1   -1    3
Please increase POLY_Dmax to at least 5 = 4 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
fi
run_test
