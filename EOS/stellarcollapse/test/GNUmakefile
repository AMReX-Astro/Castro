DEBUG := t

MKVERBOSE :=t 

COMP := gfortran

# define the location of the CASTRO root directory
CASTRO_DIR ?= ../../../
CASTRO_TOP_DIR = $(CASTRO_DIR)


# include the main Makefile stuff
include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib

# other packages needed for data_processing
Fmdirs := 


# directories containing files that are 'include'-d via Fortran
Fmincludes := 

EOS_TOP_DIR := $(CASTRO_TOP_DIR)/EOS
NETWORK_TOP_DIR := $(CASTRO_TOP_DIR)/Networks

EOS_DIR := stellarcollapse

NETWORK_DIR := general_null
GENERAL_NET_INPUTS := simple.net

MICROPHYS_CORE := $(CASTRO_TOP_DIR)/EOS

MICROPHYS_CORE += $(EOS_TOP_DIR)/$(EOS_DIR) \
                  $(NETWORK_TOP_DIR)/$(NETWORK_DIR) \
                  $(NETWORK_TOP_DIR)/

f90sources += constants_cgs.f90
f90sources += interpolate.f90
f90sources += probin.f90
f90sources += network.f90

VPATH_EXTRA_LOCATIONS += $(CASTRO_TOP_DIR)/constants
VPATH_EXTRA_LOCATIONS += $(CASTRO_TOP_DIR)/Source

Fmdirs += Util/VODE
Fmdirs += Util/BLAS


programs += testeos

all: $(pnames)

ifeq ($(findstring helmeos, $(EOS_DIR)), helmeos)
  Fmincludes := EOS/helmeos
  EOS_PATH := $(CASTRO_TOP_DIR)/EOS/$(strip $(EOS_DIR))
  all: table
endif

table:
	@if [ ! -f helm_table.dat ]; then echo Linking helm_table.dat; ln -s $(EOS_PATH)/helm_table.dat .;  fi



# any Util stuff (BLAS and VODE -- for now these seem to have a GPackage.mak)
Fmpack := $(foreach dir, $(Fmdirs), $(CASTRO_TOP_DIR)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(CASTRO_TOP_DIR)/$(dir))

# BoxLib stuff
Fmpack += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))

# Microphysics -- no GPackage.mak, so we do this manually
f90sources += $(foreach dir, $(MICROPHYS_CORE), $(notdir $(wildcard $(dir)/*.f90)))

# this puts them in the VPATH
Fmlocs += $(foreach dir, $(MICROPHYS_CORE), $(dir))

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


include $(NETWORK_TOP_DIR)/$(NETWORK_DIR)/Make.package

# probin stuff
PROBIN_PARAMETER_DIRS = ./
PROBIN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))

EXTERN_PARAMETER_DIRS += $(MICROPHYS_CORE)
EXTERN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

PROBIN_TEMPLATE := ./dummy.probin.template

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	$(BOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"


include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak

# stellarcollapse needs HDF5
# HDF5_ROOT and HDF5_INCL should be defined if not already environment variable
libraries += -lhdf5 -lhdf5_fortran -lz
LDFLAGS += -L$(HDF5_ROOT)/lib
VPATH_LOCATIONS += $(HDF5_INCL)
FINCLUDE_LOCATIONS += $(HDF5_INCL)

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
