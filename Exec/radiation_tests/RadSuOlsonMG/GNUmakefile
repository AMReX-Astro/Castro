COMP	   = gnu

DIM        = 1

USE_MPI    = TRUE
DEBUG      = FALSE

PROFILE    = FALSE

USE_RAD = TRUE


CASTRO_HOME ?= ../../..

EOS_DIR     := rad_power_law

NETWORK_DIR := general_null
NETWORK_INPUTS := gammalaw.net

OPACITY_DIR := rad_power_law

NGROUPS := 2

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro

#F90FLAGS += -x f95-cpp-input -D MG_SU_OLSON

CPPFLAGS += -DMG_SU_OLSON 
