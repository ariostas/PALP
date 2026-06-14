# Compatibility wrapper for users and tools that explicitly run
# `make -f Makefile`. The maintained make build lives in GNUmakefile.

.DEFAULT_GOAL := all

%:
	$(MAKE) -f GNUmakefile $@
