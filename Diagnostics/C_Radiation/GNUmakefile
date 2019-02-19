PRECISION = DOUBLE
PROFILE = FALSE
DEBUG = TRUE
DIM = 2

COMP = gnu

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = FALSE

USE_ACC = FALSE

ALL: gaussian_pulse.ex
# ALL: lgt_frnt1d.ex
# ALL: rad_shock.ex
# ALL: rad_source.ex
# ALL: rad_sphere.ex
# ALL: rhd_shocktube.ex

# If NEEDS_f90_SRC=TRUE, look for ${EBASE}_nd.f90
NEEDS_f90_SRC = FALSE

EOS_DIR := helmholtz

NETWORK_DIR := aprox13

Bpack   := ./Make.package
Blocs   := .
# EXTERN_SEARCH = .

CASTRO_HOME := ../..

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

include $(CASTRO_HOME)/Exec/Make.Castro

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
