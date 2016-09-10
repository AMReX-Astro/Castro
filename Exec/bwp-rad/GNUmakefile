CASTRO_HOME ?= /path/to/Castro

PRECISION        = DOUBLE
PROFILE          = FALSE
DEBUG            = FALSE
DIM              = 2
COMP	         = g++
FCOMP	         = gfortran

USE_MPI          = TRUE
USE_OMP          = FALSE

USE_GRAV         = TRUE
USE_REACT        = TRUE
USE_MODELPARSER  = TRUE

ifdef MICROPHYSICS_HOME

#This sets the EOS directory in $(MICROPHYSICS_HOME)/eos
EOS_dir     := gamma_law_general

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
Network_dir := general_null
#GENERAL_NET_INPUTS = $(MICROPHYSICS_HOME)/networks/$(Network_dir)/ignition_wdconvect.net
Network_inputs  := H_He.net

else

$(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)

endif


# power-law opacity
Opacity_dir := null

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
