SHELL=/bin/bash

UNAMEN := $(shell uname -n)

# Makefile options. We'll use the ?= method so that
# any including makefile can override these for custom tests
# (e.g. changing the number of dimensions or the EOS).

PRECISION ?= DOUBLE
PROFILE ?= FALSE

DEBUG ?= FALSE
BACKTRACE ?= FALSE
TEST ?= FALSE

USE_MPI ?= TRUE

DIM ?= 3

USE_GRAV ?= TRUE
USE_REACT ?= TRUE
USE_ROTATION ?= TRUE

DIMENSION_AGNOSTIC ?= TRUE

ifeq ($(DIM), 3)
  USE_HYBRID_MOMENTUM ?= TRUE
endif

# Set the default integrator for the burner.

INTEGRATOR_DIR ?= BS

# Default compiler.

COMP ?= gnu

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos

EOS_dir ?= helmholtz

# This sets the network directory in $(MICROPHYSICS_HOME)/networks

Network_dir ?= aprox13

# If using general_null, make sure we have a reasonable set of isotopes.

ifeq ($(Network_dir), general_null)
	GENERAL_NET_INPUTS ?= $(MICROPHYSICS_HOME)/networks/$(Network_dir)/aprox13.net
endif

# Location where our source files are coming from.

Blocs := ./
Bpack := ./Make.package

include $(CASTRO_HOME)/Exec/Make.Castro
