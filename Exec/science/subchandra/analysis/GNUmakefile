AMREX_HOME ?= ../..

DEBUG = FALSE
DIM = 2
COMP = gnu
PRECISION = DOUBLE
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE

BL_NO_FORT = TRUE

programs := fburn_weight_avg


include $(AMREX_HOME)/Tools/GNUMake/Make.defs

multiple_executables = $(addsuffix .$(machineSuffix).ex, $(programs))
default: $(multiple_executables)

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) $(multiple_executables)
