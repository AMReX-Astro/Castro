PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_ALL_CASTRO = FALSE
USE_AMR_CORE = FALSE

# define the location of the CASTRO top directory
CASTRO_HOME  ?= ../..

# This sets the EOS directory in Castro/EOS
EOS_DIR     := helmholtz

# This sets the network directory in Castro/Networks
NETWORK_DIR := general_null
NETWORK_INPUTS = ignition.net

# for the Castro runtime parameters, we don't want to use Castro::
STRUCT_USE_CASTRO_CLASS := FALSE

EXTERN_SEARCH += .

Bpack   := ./Make.package
Blocs   := .

# we explicitly want runtime_parameters.H so we have access to
# Castro's runtime parameter system
Blocs += $(CASTRO_HOME)/Source/driver $(CASTRO_HOME)/Source/hydro

include $(CASTRO_HOME)/Exec/Make.Castro
