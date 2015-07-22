
PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE


#DIM        = 1
DIM        = 2

COMP	   = gcc
FCOMP	   = gfortran

USE_MPI    = FALSE
USE_OMP    = FALSE

# define the location of the CASTRO top directory
CASTRO_DIR := ../..

ifdef MICROPHYSICS_DIR

# This sets the EOS directory in $(MICROPHYSICS_DIR)/eos
EOS_dir     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_DIR)/networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(MICROPHYSICS_DIR)/networks/$(Network_dir)/ignition.net

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_DIR appropriately)

endif

# This sets the network directory in Castro/Networks


Bpack   := ./Make.package $(CASTRO_DIR)/Util/exact_riemann/Make.package
Blocs   := . $(CASTRO_DIR)/Util/exact_riemann/

include $(CASTRO_DIR)/Exec/Make.Castro
