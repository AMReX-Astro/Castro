CASTRO_HOME := ../../..

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT  = FALSE
USE_INTERPOLATE_MODEL_VELOCITY = TRUE

USE_PROB_PARAMS = TRUE

USE_THORNADO_GPU = FALSE
USE_THORNADO_ACC = FALSE
ifeq ($(USE_THORNADO_GPU),TRUE)
  INCLUDE_LOCATIONS += $(OLCF_CUDA_ROOT)/include
  LIBRARY_LOCATIONS += $(OLCF_CUDA_ROOT)/lib64
  LIBRARIES += -lcusolver -lcublas -lcudart -lcuda
endif
ifeq ($(USE_THORNADO_ACC),TRUE)
  GENERIC_PGI_FLAGS += -acc -Minfo=accel -mcmodel=medium -ta=tesla:cc70,ptxinfo,cuda10.1
endif

USE_HDF5   = TRUE
HDF5_HOME  ?= 
ifdef OLCF_HDF5_ROOT
  HDF5_HOME := $(OLCF_HDF5_ROOT)
endif

LAPACK_HOME ?= 
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
MOMENT_CLOSURE ?= MAXIMUM_ENTROPY_CB
NEUTRINO_MATTER_SOLVER ?= FIXED_POINT_COUPLED

USE_WEAKLIB = TRUE
WEAKLIB_HOME ?= $(CASTRO_HOME)/weaklib

ifdef MICROPHYSICS_HOME

  # This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
  EOS_DIR     := weaklib

  # This sets the network directory in $(MICROPHYSICS_HOME)/networks
  NETWORK_DIR := nuclear_matter

  INTEGRATOR_DIR := VODE90
else

  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
