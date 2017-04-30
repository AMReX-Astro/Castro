AMREX_HOME ?= ../../../amrex
HERE = .

PROFILE   = FALSE

DEBUG	  = TRUE
DEBUG	  = FALSE

DIM       = 2
DIM       = 3

USE_MPI     = TRUE
USE_MPI     = FALSE

COMP      = g++

EBASE = Embiggen

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

INCLUDE_LOCATIONS  = $(AMREX_HOME)/Src/Base
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Amr
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Boundary
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata

PATHDIRS  = $(HERE)
PATHDIRS += $(AMREX_HOME)/Src/Base
PATHDIRS += $(AMREX_HOME)/Src/Amr
PATHDIRS += $(AMREX_HOME)/Src/Boundary
PATHDIRS += $(AMREX_HOME)/Src/Extern/amrdata


DEFINES += -DBL_NOLINEVALUES
DEFINES += -DBL_USEOLDREADS

CEXE_sources += $(EBASE).cpp

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package

DEFINES += -DBL_PARALLEL_IO

f90EXE_sources += probinit.f90

vpath %.f   $(PATHDIRS)
vpath %.f90 $(PATHDIRS)
vpath %.F   $(PATHDIRS)
vpath %.h   $(PATHDIRS)
vpath %.H   $(PATHDIRS)
vpath %.cpp $(PATHDIRS)

all: $(executable)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

