BOXLIB_HOME ?= ../../../BoxLib
HERE = .

PROFILE   = FALSE

DEBUG	  = TRUE
DEBUG	  = FALSE

DIM       = 2
DIM       = 3

USE_MPI     =TRUE
USE_MPI     =FALSE

COMP      = Intel
FCOMP     = Intel

EBASE = Embiggen

include $(BOXLIB_HOME)/Tools/C_mk/Make.defs

INCLUDE_LOCATIONS  = $(BOXLIB_HOME)/Src/C_BaseLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_AMRLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/C_BoundaryLib
INCLUDE_LOCATIONS += $(BOXLIB_HOME)/Src/Extern/amrdata

PATHDIRS  = $(HERE)
PATHDIRS += $(BOXLIB_HOME)/Src/C_BaseLib
PATHDIRS += $(BOXLIB_HOME)/Src/C_AMRLib
PATHDIRS += $(BOXLIB_HOME)/Src/C_BoundaryLib
PATHDIRS += $(BOXLIB_HOME)/Src/Extern/amrdata


DEFINES += -DBL_NOLINEVALUES
DEFINES += -DBL_USEOLDREADS

CEXE_sources += $(EBASE).cpp

include ./Make.package
include $(BOXLIB_HOME)/Src/C_BaseLib/Make.package
include $(BOXLIB_HOME)/Src/C_BoundaryLib/Make.package

DEFINES += -DBL_PARALLEL_IO

f90EXE_sources += probinit.f90

vpath %.f   $(PATHDIRS)
vpath %.f90 $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)

all: $(executable)

include $(BOXLIB_HOME)/Tools/C_mk/Make.rules

