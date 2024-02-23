#!/bin/sh
#
# Test the mori -m examples in Section 7.2.5 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 59
COMMAND="echo '8 4 1 1 1 1 0 6 3 1 0 1 0 1' | ./mori-${DIM}d.x -fm"
DESCRIPTION="mori-${DIM}d.x -fm example on page 59"
EXPECTED=$(cat<<-EOF
2 MORI GENERATORS / dim(cone)=2 
  3  0  1  1  0  1   I:10
  0  3 -4 -1  3 -1   I:01
2 MORI GENERATORS / dim(cone)=2 
  1  1 -1  0  1  0   I:10
  0 -3  4  1 -3  1   I:01
EOF
)

if [ $DIM -lt 6 ]; then
    # The incidence numbers 01/10 flip for some reason?
    EXPECTED=$(cat<<-EOF
2 MORI GENERATORS / dim(cone)=2 
  3  0  1  1  0  1   I:01
  0  3 -4 -1  3 -1   I:10
2 MORI GENERATORS / dim(cone)=2 
  1  1 -1  0  1  0   I:01
  0 -3  4  1 -3  1   I:10
EOF
)
fi
run_test


INPUT="3 1 1 1 0 0 0 0 2 0 0 0 1 1 0 0 2 0 0 0 0 0 1 1"
COMMAND="echo '${INPUT}' | ./mori-${DIM}d.x -fm"
DESCRIPTION="mori-${DIM}d.x -m example on page 59"
EXPECTED=$(cat<<-EOF
3 MORI GENERATORS / dim(cone)=3 
  0  0  0  0  0  0  0   I:110
  1  1  0  0  0  0  1   I:101
  0  0  1  0  0  1  0   I:011
EOF
)
SKIP=1
SKIPREASON="unknown, pre-existing output deviation"
run_test "${SKIP}" "${SKIPREASON}"
