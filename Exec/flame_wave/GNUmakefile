PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = g++
FCOMP	   = gfortran

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT  = TRUE

USE_MODELPARSER  = TRUE

USE_ROTATION  = TRUE
USE_DIFFUSION = TRUE

ifdef MICROPHYSICS_DIR

  # This sets the EOS directory in $(MICROPHYSICS_DIR)/eos
  EOS_dir     := helmholtz

  # This sets the network directory in $(MICROPHYSICS_DIR)/networks
  Network_dir := triple_alpha_plus_cago

else

  $(error Error: This problem requires the Microphysics repository. Please ensure that you have downloaded it and set $$MICROPHYSICS_DIR appropriately)

endif

Conductivity_dir := constant_opacity

Bpack   := ./Make.package
Blocs   := .

include $(CASTRO_DIR)/Exec/Make.Castro
