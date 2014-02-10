
PRECISION  = DOUBLE

DEBUG      = FALSE
#DEBUG      = TRUE

PROFILE    = FALSE

USE_MPI    = TRUE

USE_HYPRE  = TRUE
USE_RAD    = TRUE
USE_GRAV   = FALSE

# if true, sets eos to collapse -- USE_RAD must be TRUE
# if false, sets eos to gamma law
USE_NEUTRINO = FALSE
#USE_NEUTRINO = TRUE

# use radiation interpbndry
RAD_INTERP = TRUE

COMP	   = g++
FCOMP	   = gfortran
#DIM        = 1
DIM        = 2
#DIM        = 3

EOS_dir     := gamma_law_general
Network_dir := general_null
Network_inputs  := gammalaw.net
Opacity_dir := null

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro
