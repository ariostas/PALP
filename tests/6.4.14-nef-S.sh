#!/bin/sh
#
# Test the nef -S example in Section 6.4.14 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 44
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
COMMAND="echo '4 1 1 1 1' | ./nef-${DIM}d.x -f -S | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -S example on page 44"
EXPECTED=$(cat<<-EOF
4 1 1 1 1 M:35 4 N:5 4  codim=2 #part=2


#points in largest cone:
layer:  1 #p:        6 #ip:        0
layer:  2 #p:       21 #ip:        1
layer:  3 #p:       56 #ip:        6


#points in largest cone:
layer:  1 #p:       20 #ip:        0
layer:  2 #p:      105 #ip:        1
layer:  3 #p:      336 #ip:       20
H:[0] P:0 V:2 3
np=1 d:0 p:1
EOF
)
run_test
