PRECISION  = DOUBLE
PROFILE    = FALSE
TINY_PROFILE = TRUE

DEBUG      = FALSE

DIM        = 3

COMP	   = pgi

USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = TRUE

USE_GRAV   = TRUE
USE_REACT  = TRUE
USE_ROTATION  = TRUE
USE_DIFFUSION = TRUE

# define the location of the CASTRO top directory
CASTRO_HOME  := ../../../../..

GPU_COMPATIBLE_PROBLEM = TRUE

EOS_DIR     := helmholtz
NETWORK_DIR := iso7
CONDUCTIVITY_DIR := stellar

PROBLEM_DIR := ../../

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
