PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE


#DIM        = 1
DIM        = 2

COMP	   = gnu

USE_MPI    = FALSE
USE_OMP    = FALSE

USE_MHD = FALSE

USE_PROB_PARAMS = TRUE

# define the location of the CASTRO top directory
CASTRO_HOME := ../../..

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_DIR     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := general_null
NETWORK_INPUTS = ignition.net

Bpack   := ./Make.package 
#$(CASTRO_HOME)/Util/exact_riemann/Make.package
#Blocs   := . $(CASTRO_HOME)/Util/exact_riemann/
#EXTERN_CORE_ADD := $(CASTRO_HOME)/Util/exact_riemann/

include $(CASTRO_HOME)/Exec/Make.Castro
