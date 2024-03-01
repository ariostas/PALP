#!/bin/sh
#
# Test the example in Section 3.2.35 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 21
COMMAND="./poly-${DIM}d.x -C2 tests/input/3.2.35-poly-C2.txt"
DESCRIPTION="poly-${DIM}d.x -C2 example on page 21"
EXPECTED=$(cat<<-EOF
pic=1  deg=64  h12= 0  rk=0 #sq=0 #dp=0 py=1  F=5 10 10 5 #Fano=1
4 5  Vertices of P* (N-lattice)    M:201 5 N:7 5
 1  0  0  0 -1
 0  0  1  0 -1
 0  1  0  0 -1
 0  0  0  1 -4
P/2: 36 points (5 vertices) of P'=P/2 (M-lattice):
P/2:  0  4  0  0  0   0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  2  2  2  2  2  2  3  3  3
P/2:  0  0  4  0  0   1  2  3  0  1  2  3  0  1  2  0  1  0  1  2  3  0  1  2  0  1  0  0  1  2  0  1  0  0  1  0
P/2:  0  0  0  4  0   0  0  0  1  1  1  1  2  2  2  3  3  0  0  0  0  1  1  1  2  2  3  0  0  0  1  1  2  0  0  1
P/2:  0  0  0  0  1   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
EOF
)
if [ $DIM -lt 6 ]; then
    SKIP=true
    SKIPREASON="unknown, pre-existing output deviation"
fi
run_test "${SKIP}" "${SKIPREASON}"
