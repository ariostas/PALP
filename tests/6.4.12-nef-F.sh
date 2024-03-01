#!/bin/sh
#
# Test the nef -F examples in Section 6.4.12 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 41
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
INPUT="12 4 2 2 2 1 1 0 8 4 0 0 0 1 1 2"
COMMAND="echo '${INPUT}' | ./nef-${DIM}d.x -f -Lp -F3 | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -f -Lp -F3 example on page 40"
EXPECTED=$(cat<<-EOF
12 4 2 2 2 1 1 0  8 4 0 0 0 1 1 2 M:371 12 N:10 7  codim=2 #part=5
5 10  Points of Poly in N-Lattice:
    0    0    0    1    0   -1    0    0    0    0
    0    0    1    0    0   -1    0    0    0    0
   -1    4    0    0    0    0    0    1    2    0
    0   -1    0    0    1    0    0    0    0    0
   -1    2    0    0    0    1    1    1    1    0
--------------------------------------------------
    4    1    2    2    1    2    0    0    0  d=12  codim=0
    4    1    0    0    1    0    2    0    0  d=8  codim=2
    2    0    1    1    0    1    0    0    1  d=6  codim=1
    2    0    0    0    0    0    1    0    1  d=4  codim=3
    1    0    0    0    0    0    0    1    0  d=2  codim=4
--------------------------------------------- #fibrations=3
    v    v    _    _    v    _    v    p    p  cd=2  m: 35  4 n: 7 4
    v    _    v    v    _    v    v    p    v  cd=1  m:117  9 n: 8 6
    v    _    _    _    _    _    v    p    v  cd=3  m:  9  3 n: 5 3
H:4 58 [-108] P:1 V:0 2   (6 6) (4 4) (3 3) (2 2) (1 1)
H:3 65 [-124] P:2 V:0 2 3   (8 4) (4 4) (4 2) (2 2) (1 1)
H:3 83 [-160] P:3 V:3 5   (4 8) (0 8) (2 4) (0 4) (0 2)
np=3 d:0 p:2
EOF
)
if [ $DIM -lt 6 ]; then
   EXPECTED="Please increase POLY_Dmax to at least 6 = 5 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
fi
run_test
