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
# in to GNU Make. There's a special case for an empty DIMENSION as an
# indicator that we should not override the value of POLY_Dmax in
# Global.h. This is used to build the foo.x programs using the value
# in Global.h, since editing Global.h has long been documented as the
# way to change POLY_Dmax in foo.x.
ifeq ($(2),)
%-$(2)d.o: %.c
	$(COMPILE.c) -o $$@ $$<
else
%-$(2)d.o: %.c
	$(COMPILE.c) -DPOLY_Dmax=$(2) -o $$@ $$<
endif

# Link the program foo-Nd.x from foo-Nd.o, OBJECTS_N, and foo_OBJ_N.
# The LINK.c macro is built in to GNU Make.
$(1)-$(2)d.x: $(1)-$(2)d.o $$(OBJECTS_$(2)) $$($(1)_OBJ_$(2))
	$(LINK.c) -o $$@ $$^

# Add foo-Nd.x to the "all-dims" target
all-dims: $(1)-$(2)d.x

# Specify some additional dependencies (beyond the corresponding *.c file)
# for our *.o files.
Coord-$(2)d.o:         Rat.h Global.h
Polynf-$(2)d.o:        Rat.h Global.h
Rat-$(2)d.o:           Rat.h Global.h
Subpoly-$(2)d.o:       Rat.h Subpoly.h Global.h
Subadd-$(2)d.o:        Subpoly.h Global.h
Vertex-$(2)d.o:        Rat.h Global.h
Subdb-$(2)d.o:         Subpoly.h Global.h
LG-$(2)d.o:            Rat.h LG.h Global.h

E_Poly-$(2)d.o:        Nef.h Rat.h Global.h
Nefpart-$(2)d.o:       Nef.h Global.h

MoriCone-$(2)d.o:      Rat.h Mori.h Global.h
SingularInput-$(2)d.o: Mori.h Global.h

poly-$(2)d.o:          LG.h Global.h
class-$(2)d.o:         Subpoly.h Global.h
cws-$(2)d.o:           LG.h Rat.h Global.h
nef-$(2)d.o:           Nef.h LG.h Global.h
mori-$(2)d.o:          LG.h Mori.h Global.h
endef

# Call the PROG_DIM_template once for each PROGRAM "p" and
# DIMENSION "d".
$(foreach p,$(PROGRAMS),$(foreach d,$(DIMENSIONS),\
  $(eval $(call PROG_DIM_template,$(p),$(d)))\
))

# Typically (i.e. by default) foo.x will be identical to foo-6d.x,
# since Global.h defines POLY_Dmax=6. However, the procedure to change
# POLY_Dmax has long been documented as "edit Global.h and re-run
# make." If we simply copy foo-6d.x to foo.x, then changes to Global.h
# will have no effect on the resulting executables. Instead, we build
# a separate copy of each program with DIMENSION set to the empty
# string. This alerts the build rule to omit the -DPOLY_Dmax line that
# overrides the value in Global.h. The resulting foo-d.x executables
# can then be copied to foo.x so that the documented procedure still
# works. The special INTERMEDIATE rule avoids rebuilding foo-d.x after
# we've moved it to foo.x.
$(foreach p,$(PROGRAMS),\
  $(eval $(call PROG_DIM_template,$(p),))\
)
.INTERMEDIATE: $(foreach p,$(PROGRAMS),$(p)-d.x)
%.x: %-d.x
	mv $< $@


# For lack of anything less silly, define a newline this way.
# We need it to run multiple commands in a $(foreach) loop.
# The $(blank) is apparently required to keep it from eating
# the newline.
blank :=
define newline

$(blank)
endef

# Testing
# #######
#
# Each test script should take the dimension DIM from the environment
# with default 6. Optionally, the LONG variable can be set to a
# non-null value to enable the tests that take a long time. Some
# documentation on this test "protocol" is contained in tests/README.
#
TESTS = $(wildcard tests/*.sh)

# The test suite depends on all of the dimension-optimized programs.
# The main "make check" routine runs each of the tests/*.sh scripts
# with each applicable value of DIM.
.PHONY: check
check: $(foreach p,$(PROGRAMS),$(foreach d,$(DIMENSIONS),$(p)-$(d)d.x))
	$(foreach t,$(TESTS),$(foreach d,$(DIMENSIONS),\
	$(newline)@DIM=$(d) $(t)\
	))

.PHONY: checklong
checklong:
	LONG=1 $(MAKE) check
