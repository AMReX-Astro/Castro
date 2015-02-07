COMP	   = g++
FCOMP	   = gfortran

DIM        = 1

DEBUG      = TRUE
DEBUG      = FALSE

PROFILE    = FALSE

USE_NEUTRINO = FALSE

EOS_dir     := gamma_law_general
Network_dir := general_null
Network_inputs  := gammalaw.net
Opacity_dir := null

Bpack   := ./Make.package
Blocs   := .

include ../Make.Castro

#F90FLAGS += -x f95-cpp-input -D MG_SU_OLSON

CPPFLAGS += -DMG_SU_OLSON 
