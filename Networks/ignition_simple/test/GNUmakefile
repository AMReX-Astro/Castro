NDEBUG :=  
NDEBUG := 
MPI    :=
OMP    :=

COMP := Intel

MKVERBOSE := t

FPARALLEL := ../../../..

include $(FPARALLEL)/mk/GMakedefs.mak

programs += testburn
programs += eval
programs += testjacobian

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/extern/helmeos/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/helmeos
FINCLUDE_LOCATIONS += $(FPARALLEL)/extern/helmeos

include $(FPARALLEL)/extern/networks/ignition_simple/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/networks/ignition_simple/

include $(FPARALLEL)/extern/VODE/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/VODE

all: $(pnames)

include $(FPARALLEL)/mk/GMakerules.mak

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif



