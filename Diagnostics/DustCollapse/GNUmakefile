BOXLIB_HOME ?= /path/to/BoxLib

NDEBUG := t
MPI    :=
OMP    :=
MKVERBOSE :=t 
COMP   := gfortran

# include the main Makefile stuff
include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# define the packages to build these tools
Fmdirs := Src/F_BaseLib 

Fmpack := $(foreach dir, $(Fmdirs),     $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs),     $(BOXLIB_HOME)/$(dir))
Fmincs := $(foreach dir, $(Fmdirs),     $(BOXLIB_HOME)/$(dir))

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)

programs += fdustcollapse1d
programs += fdustcollapse2d
programs += fdustcollapse3d

all: $(pnames)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else	
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif
