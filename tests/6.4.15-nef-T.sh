#!/bin/sh
#
# Test the nef -S -T example in Section 6.4.15 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 44-45
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
COMMAND="echo '4 1 1 1 1' | ./nef-${DIM}d.x -f -S -T | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -S -T example on pages 44-45"
EXPECTED=$(cat<<-EOF
4 1 1 1 1 M:35 4 N:5 4  codim=2 #part=2


#points in largest cone:
layer:  1 #p:        6 #ip:        0
layer:  2 #p:       21 #ip:        1
layer:  3 #p:       56 #ip:        6
layer:  4 #p:      125 #ip:       21
layer:  5 #p:      246 #ip:       56


#points in largest cone:
layer:  1 #p:       20 #ip:        0
layer:  2 #p:      105 #ip:        1
layer:  3 #p:      336 #ip:       20
layer:  4 #p:      825 #ip:      105
layer:  5 #p:     1716 #ip:      336
H:[0] P:0 V:2 3
np=1 d:0 p:1
EOF
)
run_test
