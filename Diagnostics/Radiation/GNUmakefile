PRECISION = DOUBLE
PROFILE = FALSE
DEBUG = TRUE
DIM = 2

COMP = gnu

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = FALSE

ALL: radhelp

radhelp:
	$(error You need to specify a target, e.g., 'make DIM=1 rad_sphere.ex')

# If NEEDS_f90_SRC=TRUE, look for ${EBASE}_nd.f90
NEEDS_f90_SRC = FALSE

EOS_DIR := helmholtz

NETWORK_DIR := aprox13

Bpack   := ./Make.package
Blocs   := .
# EXTERN_SEARCH = .

CASTRO_HOME ?= ../..

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

include $(CASTRO_HOME)/Exec/Make.Castro

ifeq ($(MAKECMDGOALS),rad_sphere.ex)
  CEXE_sources += rad_sphere.cpp
endif

ifeq ($(MAKECMDGOALS),gaussian_pulse.ex)
  CEXE_sources += gaussian_pulse.cpp
endif

ifeq ($(MAKECMDGOALS),lgt_frnt1d.ex)
  CEXE_sources += lgt_frnt1d.cpp
endif

ifeq ($(MAKECMDGOALS),rad_shock.ex)
  CEXE_sources += rad_shock.cpp
endif

ifeq ($(MAKECMDGOALS),rad_source.ex)
  CEXE_sources += rad_source.cpp
endif

ifeq ($(MAKECMDGOALS),rhd_shocktube.ex)
  CEXE_sources += rhd_shocktube.cpp
endif

gaussian_pulse.ex: $(objForExecs)
	@echo Linking $@ ...
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)

lgt_frnt1d.ex: $(objForExecs)
	@echo Linking $@ ...
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)

rad_shock.ex: $(objForExecs)
	@echo Linking $@ ...
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)

rad_source.ex: $(objForExecs)
	@echo Linking $@ ...
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)

rad_sphere.ex: $(objForExecs)
	@echo Linking $@ ...
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)

rhd_shocktube.ex: $(objForExecs)
	@echo Linking $@ ...
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)
