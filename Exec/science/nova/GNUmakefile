PRECISION        = DOUBLE
PROFILE          = FALSE
DEBUG            = FALSE
DIM              = 2

COMP	         = gnu

USE_MPI          = TRUE
USE_OMP          = FALSE
USE_CUDA         = TRUE

USE_REACT        = TRUE
USE_DIFFUSION    = TRUE
USE_GRAV         = TRUE

CASTRO_HOME := ../../..

USE_CXX_MODEL_PARSER = TRUE

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := nova2

# Thi sets the conductivity EOS directory in $(MICROPHYSICS_HOME)/conductivity
CONDUCTIVITY_DIR := stellar

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_HOME)/Exec/Make.Castro
