PRECISION ?= DOUBLE
PROFILE ?= FALSE

DEBUG ?= FALSE
BACKTRACE ?= FALSE
TEST ?= FALSE

USE_MPI ?= TRUE
USE_OMP ?= FALSE

DIM ?= 3

USE_GRAV ?= TRUE
USE_REACT ?= TRUE
USE_ROTATION ?= TRUE

# We want to use the dimension agnostic version of certain source files.

GPU_COMPATIBLE_PROBLEM = TRUE

# Turn on hybrid momentum, but note that it only is used in 3D at present.

ifeq ($(DIM), 3)
  USE_HYBRID_MOMENTUM ?= TRUE
endif

# Set the default integrator for the burner.

INTEGRATOR_DIR ?= VODE90

# Default compiler.

COMP ?= gnu

ifndef MICROPHYSICS_HOME
  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)
endif

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos

EOS_DIR ?= helmholtz

# This sets the network directory in $(MICROPHYSICS_HOME)/networks

NETWORK_DIR ?= aprox13

# If using general_null, make sure we have a reasonable set of isotopes.

ifeq ($(NETWORK_DIR), general_null)
	GENERAL_NET_INPUTS = $(MICROPHYSICS_HOME)/networks/$(NETWORK_DIR)/aprox13.net
endif

# Location where the source files are coming from.

Blocs := $(CASTRO_HOME)/Exec/science/wdmerger/
Bpack := $(CASTRO_HOME)/Exec/science/wdmerger/Make.package

include $(CASTRO_HOME)/Exec/Make.Castro
