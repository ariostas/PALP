#!/bin/sh
#
# Test the mori -t example in Section 7.2.11 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Page 63
COMMAND="echo '8 4 1 1 1 1 0 6 3 1 0 1 0 1' | ./mori-${DIM}d.x -ft | head -n37"
DESCRIPTION="mori-${DIM}d.x -ft example on page 63"
EXPECTED=$(cat<<-EOF
SINGULAR -> triple intersection numbers:
d6^3->2,
d5*d6^2->2,
d4*d6^2->2,
d3*d6^2->0,
d2*d6^2->2,
d1*d6^2->8,
d5^2*d6->0,
d4*d5*d6->2,
d3*d5*d6->2,
d2*d5*d6->0,
d1*d5*d6->6,
d4^2*d6->2,
d3*d4*d6->0,
d2*d4*d6->2,
d1*d4*d6->8,
d3^2*d6->-2,
d2*d3*d6->2,
d1*d3*d6->2,
d2^2*d6->0,
d1*d2*d6->6,
d1^2*d6->30,
d5^3->0,
d4*d5^2->0,
d3*d5^2->0,
d2*d5^2->0,
d1*d5^2->0,
d4^2*d5->2,
d3*d4*d5->2,
d2*d4*d5->0,
d1*d4*d5->6,
d3^2*d5->2,
d2*d3*d5->0,
d1*d3*d5->6,
d2^2*d5->0,
d1*d2*d5->0,
d1^2*d5->18,
EOF
)
run_test
