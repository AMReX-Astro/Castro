PRECISION = DOUBLE
PROFILE = FALSE
DEBUG = TRUE
DIM = 3

COMP = gnu

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = FALSE

# programs to be compiled
ALL: sedov_$(DIM)d.ex

EOS_DIR := helmholtz

NETWORK_DIR := aprox13

Bpack   := ./Make.package
Blocs   := .
# EXTERN_SEARCH = .

CASTRO_HOME ?= ../..

#INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
#include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
#vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

include $(CASTRO_HOME)/Exec/Make.Castro

sedov_$(DIM)d.ex: $(objForExecs)
	@echo Linking $@ ...
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)
