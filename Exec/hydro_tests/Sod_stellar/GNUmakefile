PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE


#DIM        = 1
DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE

USE_MHD = FALSE


# define the location of the CASTRO top directory
CASTRO_HOME ?= ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := helmholtz

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := general_null
NETWORK_INPUTS = ignition.net

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
#$(CASTRO_HOME)/Util/exact_riemann/Make.package
#Blocs   := $(PROBLEM_DIR) $(CASTRO_HOME)/Util/exact_riemann/
#EXTERN_CORE_ADD := $(CASTRO_HOME)/Util/exact_riemann/

include $(CASTRO_HOME)/Exec/Make.Castro
