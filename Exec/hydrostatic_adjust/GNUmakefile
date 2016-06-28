PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 1

COMP	   = Intel
FCOMP	   = Intel

COMP	   = g++
FCOMP	   = gfortran

USE_MPI    = TRUE

USE_GRAV   = TRUE

CASTRO_HOME = ../..

USE_MODELPARSER = TRUE

ifdef MICROPHYSICS_HOME

# This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_dir     := helmholtz

# This sets the EOS directory in $(MICROPHYSICS_HOME)/networks
Network_dir := general_null
GENERAL_NET_INPUTS = $(MICROPHYSICS_HOME)/networks/$(Network_dir)/ignition.net

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif

Bpack	:= ./Make.package
Blocs	:= .

include $(CASTRO_HOME)/Exec/Make.Castro
