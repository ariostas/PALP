#!/bin/sh
#
# Test the nef -H examples in Section 6.4.4 of the PALP manual
#
. tests/lib/run-test.sh

: ${DIM:=6}

# Pages 35
# The "sed" command strips the unpredictable timing information
# from the ends of the lines.
COMMAND="echo '7 1 1 1 1 1 1 1' | ./nef-${DIM}d.x -f -H | sed 's/ *[0-9]*sec.*//g'"
DESCRIPTION="nef-${DIM}d.x -H example on page 35"
EXPECTED=$(cat<<-EOF
7 1 1 1 1 1 1 1 M:1716 7 N:8 7  codim=2 #part=3


                                    h 0 0   

                               h 1 0      h 0 1   

                          h 2 0      h 1 1      h 0 2   

                     h 3 0      h 2 1      h 1 2      h 0 3   

                h 4 0      h 3 1      h 2 2      h 1 3      h 0 4   

                     h 4 1      h 3 2      h 2 3      h 1 4   

                          h 4 2      h 3 3      h 2 4   

                               h 4 3      h 3 4   

                                    h 4 4   



                                       1

                                  0         0

                             0         1         0

                        0         0         0         0

                   1       237       996       237         1

                        0         0         0         0

                             0         1         0

                                  0         0

                                       1




                                    h 0 0   

                               h 1 0      h 0 1   

                          h 2 0      h 1 1      h 0 2   

                     h 3 0      h 2 1      h 1 2      h 0 3   

                h 4 0      h 3 1      h 2 2      h 1 3      h 0 4   

                     h 4 1      h 3 2      h 2 3      h 1 4   

                          h 4 2      h 3 3      h 2 4   

                               h 4 3      h 3 4   

                                    h 4 4   



                                       1

                                  0         0

                             0         1         0

                        0         0         0         0

                   1       356      1472       356         1

                        0         0         0         0

                             0         1         0

                                  0         0

                                       1


np=2 d:0 p:1
EOF
)
if [ $DIM -lt 7 ]; then
   EXPECTED="Please increase POLY_Dmax to at least 7 = 6 + 2 - 1
(POLY_Dmax >= dim N + codim - 1 is required)"
else
    if [ -z "${LONG}" ]; then
	# This test takes forever to run
	SKIP=true
	SKIPREASON="long-running test requires make checklong"
    fi
fi
run_test "${SKIP}" "${SKIPREASON}"
