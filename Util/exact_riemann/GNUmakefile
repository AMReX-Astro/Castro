DEBUG := t

MKVERBOSE :=t 

COMP := gfortran

# define the location of the CASTRO root directory
CASTRO_DIR ?= ../../
CASTRO_TOP_DIR = $(CASTRO_DIR)


# include the main Makefile stuff
include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib

# other packages needed for data_processing
Fmdirs := 


# directories containing files that are 'include'-d via Fortran
Fmincludes :=

ifndef MICROPHYSICS_DIR
$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_DIR appropriately)
endif

USE_REACT := FALSE
EOS_TOP_DIR := $(MICROPHYSICS_DIR)/eos
NETWORK_TOP_DIR := $(MICROPHYSICS_DIR)/networks

EOS_DIR := gamma_law_general
#EOS_DIR := helmholtz

NETWORK_DIR := ignition_simple

MICROPHYS_CORE := eos networks
MICROPHYS_CORE += eos/$(EOS_DIR) networks/$(NETWORK_DIR)

f90sources += probin.f90

Fmdirs += Util/VODE
Fmdirs += Util/BLAS


programs += exact_riemann

all: $(pnames)

ifeq ($(findstring helmholtz, $(EOS_DIR)), helmholtz)
  all: helm_table
endif

helm_table:
	@if [ ! -f helm_table.dat ]; then echo Linking helm_table.dat; ln -s $(EOS_TOP_DIR)/$(EOS_DIR)/helm_table.dat .;  fi



# any Util stuff (BLAS and VODE -- for now these seem to have a GPackage.mak)
Fmpack := $(foreach dir, $(Fmdirs), $(CASTRO_TOP_DIR)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(CASTRO_TOP_DIR)/$(dir))

# BoxLib stuff
Fmpack += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))

# EOS/Network
Fmpack += $(foreach dir, $(MICROPHYS_CORE), $(MICROPHYSICS_DIR)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(MICROPHYS_CORE), $(MICROPHYSICS_DIR)/$(dir))

# CASTRO-specific stuff that the Microphysics directories rely on
Fmlocs += $(CASTRO_TOP_DIR)/EOS
f90sources += eos.f90
f90sources += eos_data.f90
f90sources += eos_type.f90

Fmlocs += $(CASTRO_TOP_DIR)/Networks
f90sources += network.f90

Fmlocs += $(CASTRO_TOP_DIR)/constants
f90sources += constants_cgs.f90

Fmlocs += $(CASTRO_TOP_DIR)/Source/Src_nd/
f90sources += meth_params.f90

# any include files
Fmincs := $(foreach dir, $(Fmincludes), $(CASTRO_TOP_DIR)/$(dir))

# current directory
Fmpack += ./GPackage.mak

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs) $(VPATH_EXTRA_LOCATIONS)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)



# probin stuff
PROBIN_PARAMETER_DIRS = ./
PROBIN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))

EXTERN_PARAMETER_DIRS += $(foreach dir, $(MICROPHYS_CORE), $(MICROPHYSICS_DIR)/$(dir))
EXTERN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

PROBIN_TEMPLATE := ./dummy.probin.template

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	$(BOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"


include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak


%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else	
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif


# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will print out the
# value.
print-%: ; @echo $* is $($*)



clean::
	$(RM) probin.f90
