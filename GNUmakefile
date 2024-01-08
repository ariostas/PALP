# 				M A K E  F I L E

#   main programs:	 class.c  cws.c  poly.c  nef.c  mori.c

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

# The list of programs to build, without the ".x" extension. Having
# these in a variable makes it easy to loop through them. Likewise,
# the list of dimensions (POLY_Dmax) that we want to build executables
# for.
PROGRAMS = poly class cws nef mori
DIMENSIONS = 4 5 6 11

# The "all" target builds only the 6d-optimized versions that have
# historically been built by running "make".
.PHONY: all
all: $(foreach p,$(PROGRAMS),$(p).x)

# The "all-dims" target, however, builds each PROGRAM, once, optimized
# for each dimension listed in DIMENSIONS. Here we have "all-dims"
# depend only on the traditional foo.x names, but the template below
# will add foo-4d.x, foo-5d.x, etc. to the list of prerequisites.
.PHONY: all-dims
all-dims: all


.PHONY: clean
clean:	;	rm -f *.o

.PHONY: cleanall
cleanall: ;	rm -f *.o *.x palp_* core


# Build foo.x by copying foo-6d.x to foo.x. This wastes a little bit
# of space, but avoids any obscure problems that might arise from
# using links instead
%.x: %-6d.x
	cp $< $@

define PROG_DIM_template =
#
# Define separate build rules for every combination of PROGRAMS and
# DIMENSIONS. This really is necessary: we can't reuse an object file
# that was compiled with (say) POLY_Dmax=4 to link the executable
# foo-11d.x, because then foo-11d.x will just wind up with the code
# for dimension <= 4. And that's the best case: mixing and matching
# POLY_Dmax across multiple files could easily cause a crash.
#
# Arguments:
#
#   $(1) - program name, e.g. "poly" or "mori"
#   $(2) - the current value of POLY_Dmax, e.g. "4" or "11"
#

# A list of common objects needed by all executables of this dimension
OBJECTS_$(2)    = Coord-$(2)d.o Rat-$(2)d.o Vertex-$(2)d.o Polynf-$(2)d.o

# List the additional objects needed by the individual programs of
# this dimension
poly_OBJ_$(2)   = LG-$(2)d.o
class_OBJ_$(2)  = Subpoly-$(2)d.o Subadd-$(2)d.o Subdb-$(2)d.o
cws_OBJ_$(2)    = LG-$(2)d.o
nef_OBJ_$(2)    = E_Poly-$(2)d.o Nefpart-$(2)d.o LG-$(2)d.o
mori_OBJ_$(2)   = MoriCone-$(2)d.o SingularInput-$(2)d.o LG-$(2)d.o

# Build the object foo-Nd.o from foo.c. The COMPILE.c macro is built
# in to GNU Make.
%-$(2)d.o: %.c
	$(COMPILE.c) -DPOLY_Dmax=$(2) -o $$@ $$<

# Link the program foo-Nd.x from foo-Nd.o, OBJECTS_N, and foo_OBJ_N.
# The LINK.c macro is built in to GNU Make.
$(1)-$(2)d.x: $(1)-$(2)d.o $$(OBJECTS_$(2)) $$($(1)_OBJ_$(2))
	$(LINK.c) -o $$@ $$^

# Add foo-Nd.x to the "all-dims" target
all-dims: $(1)-$(2)d.x

# Specify some additional dependencies (beyond the corresponding *.c file)
# for our *.o files.
Coord-$(2)d.o:         Rat.h
Polynf-$(2)d.o:        Rat.h
Rat-$(2)d.o:           Rat.h
Subpoly-$(2)d.o:       Rat.h Subpoly.h
Subadd-$(2)d.o:        Subpoly.h
Vertex-$(2)d.o:        Rat.h
Subdb-$(2)d.o:         Subpoly.h
LG-$(2)d.o:            Rat.h LG.h

E_Poly-$(2)d.o:        Nef.h Rat.h
Nefpart-$(2)d.o:       Nef.h

MoriCone-$(2)d.o:      Rat.h Mori.h
SingularInput-$(2)d.o: Mori.h

poly-$(2)d.o:          LG.h
class-$(2)d.o:         Subpoly.h
cws-$(2)d.o:           LG.h Rat.h
nef-$(2)d.o:           Nef.h LG.h
mori-$(2)d.o:          LG.h Mori.h
endef

# All object files should be rebuilt if Global.h changes
%.o: Global.h

# Call the PROG_DIM_template once for each PROGRAM "p" and
# DIMENSION "d".
$(foreach p,$(PROGRAMS),$(foreach d,$(DIMENSIONS),\
  $(eval $(call PROG_DIM_template,$(p),$(d)))\
))
