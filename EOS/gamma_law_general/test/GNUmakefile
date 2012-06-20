NDEBUG :=  
NDEBUG := t
MPI    :=
OMP    :=

COMP := Intel

MKVERBOSE := t

FPARALLEL := ../../../..

include $(FPARALLEL)/mk/GMakedefs.mak

f90sources += testeos.f90

# define the packages to build this problem
#   note: we only need to include the problem's directory if there are 
#   unique files there (as specified in a GPackage.mak).  The problem 
#   directory is always placed at the start of the vpath by the 
#   GMakerules.mak.
Fmdirs := boxlib \
          extern/EOS/gamma_law_general \
          extern/networks/ignition \
          extern/VODE \
          extern/constants

# directories containing files that are 'include'-d via Fortran
Fmincludes := extern/helmeos


Fmpack := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir))
Fmincs := $(foreach dir, $(Fmincludes), $(FPARALLEL)/$(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


# define the build instructions for the executable
testeos.$(suf).exe: $(objects)
	$(LINK.f90) -o testeos.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

