# 				M A K E  F I L E

#   main programs:	 class.c  cws.c  poly.c  nef.c  mori.c

SOURCES= Coord.c Rat.c Vertex.c Polynf.c
OBJECTS= $(SOURCES:.c=.o)

CLASS_SRC= Subpoly.c Subadd.c Subdb.c
CLASS_OBJ= $(CLASS_SRC:.c=.o)

NEF_SRC= E_Poly.c Nefpart.c LG.c
NEF_OBJ= $(NEF_SRC:.c=.o)

MORI_SRC= MoriCone.c SingularInput.c
MORI_OBJ= $(MORI_SRC:.c=.o)


CC ?= gcc

CPPFLAGS += -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
CFLAGS ?= -O3 -g -W -Wall
# CFLAGS=-O3 -g				      # add -g for GNU debugger gdb
# CFLAGS=-Ofast -O3 -mips4 -n32		      # SGI / 32 bit
# CFLAGS=-Ofast -O3 -mips4 -64                # SGI / 64 bit
# CFLAGS=-O3 -fast -non_shared -arch pca56    # statically linked for alpha_PC

#   targets : dependencies ; command
#             command
#             ...

all:	poly class cws nef mori

clean:	;	rm -f *.o

cleanall: ;	rm -f *.o *.x palp_* core


poly:	poly.o $(OBJECTS) LG.o Global.h LG.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o poly.x poly.o $(OBJECTS) LG.o

class:	class.o $(OBJECTS) $(CLASS_OBJ) Global.h Subpoly.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o class.x \
                class.o $(OBJECTS) $(CLASS_OBJ)

cws:    cws.o $(OBJECTS) LG.o Global.h LG.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o cws.x \
                cws.o $(OBJECTS) LG.o

nef:    nef.o $(OBJECTS) $(NEF_OBJ) Global.h 
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o nef.x \
                nef.o $(OBJECTS) $(NEF_OBJ)

mori:   mori.o  $(OBJECTS) $(MORI_OBJ) LG.o Mori.h 
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o mori.x \
                mori.o $(OBJECTS) $(MORI_OBJ) LG.o



#			     D E P E N D E N C I E S

Coord.o:        Global.h Rat.h 
Polynf.o:		Global.h Rat.h
Rat.o:          Global.h Rat.h 
Subpoly.o:      Global.h Rat.h Subpoly.h  
Subadd.o:		Global.h Subpoly.h
Vertex.o:       Global.h Rat.h  
Subdb.o:		Global.h Subpoly.h
LG.o:           Global.h Rat.h LG.h

E_Poly.o:       Global.h Nef.h Rat.h
Nefpart.o:		Global.h Nef.h

MoriCone.o:      Global.h Rat.h Mori.h
SingularInput.o: Global.h Mori.h

poly.o:         Global.h LG.h
class.o:		Global.h Subpoly.h
cws.o:			Global.h LG.h Rat.h
nef.o:          Global.h Nef.h LG.h
mori.o:     	Global.h LG.h Mori.h
