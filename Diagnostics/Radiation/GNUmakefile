NDEBUG := t
MPI    :=
OMP    :=
MKVERBOSE :=t 

COMP := gfortran


# include the main Makefile stuff
include $(AMREX_HOME)/Tools/F_mk/GMakedefs.mak

# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib

# BoxLib stuff
Fmpack += $(foreach dir, $(BOXLIB_CORE), $(AMREX_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(AMREX_HOME)/$(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)



programs += flgt_frnt1d.f90


all: $(pnames)

include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else	
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif



