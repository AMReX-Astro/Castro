DEBUG      = FALSE

DIM        = 3

COMP       = gnu

USE_MPI    = TRUE

USE_RAD    = TRUE

CASTRO_HOME ?= ../../..

EOS_DIR     := gamma_law

NETWORK_DIR := general_null
NETWORK_INPUTS := gammalaw.net

Opacity_dir := rad_power_law

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
