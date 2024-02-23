#!/bin/sh
#
# Test the nef -s -Lv example in Section 6.4.16 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 45
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
INPUT="3 1 1 1 0 0 0 3 0 0 0 1 1 1"
COMMAND="echo '${INPUT}' | ./nef-${DIM}d.x -f -s -Lv | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -s -Lv example on page 45"
EXPECTED=$(cat<<-EOF
3 1 1 1 0 0 0  3 0 0 0 1 1 1 M:100 9 N:7 6  codim=2 #part=31
4 6 Vertices in N-lattice:
    0    0    0    1    0   -1
    0    0    1    0    0   -1
   -1    0    0    0    1    0
   -1    1    0    0    0    0
------------------------------
    1    1    0    0    1    0  d=3  codim=2
    0    0    1    1    0    1  d=3  codim=2
H:20 [24] P:2 V:4 5   (1 2) (1 2)
H:20 [24] P:4 V:0 5   (1 2) (1 2)
H:20 [24] P:5 V:0 4   (2 1) (0 3)
H:20 [24] P:6 V:0 4 5   (2 1) (1 2)
H:20 [24] P:8 V:1 5   (1 2) (1 2)
H:20 [24] P:9 V:1 4   (2 1) (0 3)
H:20 [24] P:10 V:1 4 5   (2 1) (1 2)
H:20 [24] P:11 V:0 1   (2 1) (0 3)
H:20 [24] P:12 V:0 1 5   (2 1) (1 2)
H:20 [24] P:14 V:2 3   (0 3) (2 1)
H:20 [24] P:16 V:2 5   (0 3) (2 1)
H:20 [24] P:17 V:2 4   (1 2) (1 2)
H:20 [24] P:18 V:2 4 5   (1 2) (2 1)
H:20 [24] P:19 V:0 2   (1 2) (1 2)
H:20 [24] P:20 V:0 2 5   (1 2) (2 1)
H:20 [24] P:21 V:0 2 4   (2 1) (1 2)
H:20 [24] P:22 V:1 3   (1 2) (1 2)
H:20 [24] P:23 V:1 2   (1 2) (1 2)
H:20 [24] P:24 V:1 2 5   (1 2) (2 1)
H:20 [24] P:25 V:1 2 4   (2 1) (1 2)
H:20 [24] P:26 V:0 3   (1 2) (1 2)
H:20 [24] P:27 V:0 1 2   (2 1) (1 2)
H:20 [24] P:28 V:3 4   (1 2) (1 2)
H:20 [24] P:29 V:3 5   (0 3) (2 1)
np=24 d:1 p:6
EOF
)
if [ $DIM -lt 5 ]; then
   EXPECTED="Please increase POLY_Dmax to at least 5 = 4 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
elif [ $DIM -eq 5 ]; then
    SKIP=1
    SKIPREASON="Unknown, pre-existing output deviation"
fi
run_test "${SKIP}" "${SKIPREASON}"
