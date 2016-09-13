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
CASTRO_HOME := ../../..

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_dir     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_HOME)/networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(MICROPHYSICS_HOME)/networks/$(Network_dir)/ignition.net

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

# This sets the network directory in Castro/Networks


Bpack   := ./Make.package 
#$(CASTRO_HOME)/Util/exact_riemann/Make.package
#Blocs   := . $(CASTRO_HOME)/Util/exact_riemann/
#EXTERN_CORE_ADD := $(CASTRO_HOME)/Util/exact_riemann/

include $(CASTRO_HOME)/Exec/Make.Castro
