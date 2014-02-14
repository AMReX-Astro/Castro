NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

FPARALLEL := ../../..

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

f90sources += testeos.f90

EOS_DIR := EOS/helmeos

NETWORK_DIR := Networks/general_null
GENERAL_NET_INPUTS := $(FPARALLEL)/$(NETWORK_DIR)/ignition.net

BOXLIB_CORE := Src/F_BaseLib

ALL: testeos.$(suf).exe

# define the packages to build this problem
#   note: we only need to include the problem's directory if there are 
#   unique files there (as specified in a GPackage.mak).  The problem 
#   directory is always placed at the start of the vpath by the 
#   GMakerules.mak.
Fmdirs := $(EOS_DIR) \
	  $(NETWORK_DIR)

# directories containing files that are 'include'-d via Fortran
#Fmincludes := extern/helmeos


Fmpack := $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))
Fmincs :=

Fmpack += $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir)/Make.package)
Fmlocs += $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir))
Fmincs += $(foreach dir, $(Fmincludes), $(FPARALLEL)/$(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


# define the build instructions for the executable
testeos.$(suf).exe: $(objects)
	$(LINK.f90) -o testeos.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak

