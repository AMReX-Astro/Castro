PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3
DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE
USE_GPU    = FALSE
USE_ACC    = FALSE

USE_HDF5   = TRUE
ifdef OLCF_HDF5_ROOT
  HDF5_HOME := $(OLCF_HDF5_ROOT)
endif
ifeq ($(HOST),manda)
  HDF5_HOME := /home/almgren/bin/miniconda2
endif

ifdef OLCF_NETLIB_LAPACK_ROOT
  LAPACK_HOME := $(OLCF_NETLIB_LAPACK_ROOT)/lib64
endif
ifdef OLCF_ESSL_ROOT
  INCLUDE_LOCATIONS += $(OLCF_ESSL_ROOT)/include
  LIBRARY_LOCATIONS += $(OLCF_ESSL_ROOT)/lib64
  LIBRARIES += -lessl
endif

USE_THORNADO = TRUE
THORNADO_HOME ?= $(CASTRO_HOME)/thornado
MOMENT_CLOSURE ?= MINERBO
NEUTRINO_MATTER_SOLVER ?= EMAB

USE_WEAKLIB = FALSE
WEAKLIB_HOME ?= $(CASTRO_HOME)/weaklib

# define the location of the CASTRO top directory
CASTRO_HOME ?= ../../..
MICROPHYSICS_HOME ?= ../../../../Microphysics

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := general_null

GENERAL_NET_INPUTS = ./ignition.net

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

# This sets the network directory in Castro/Networks

Bpack   := ./Make.package 

include $(CASTRO_HOME)/Exec/Make.Castro
