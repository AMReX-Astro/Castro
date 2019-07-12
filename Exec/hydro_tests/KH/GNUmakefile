# Define the location of the CASTRO top directory,
# if not already defined by an environment variable.

CASTRO_HOME := ../../..

# Location of this directory. Useful if
# you're trying to compile this from another location.

TEST_DIR = $(CASTRO_HOME)/Exec/hydro_tests/KH

PRECISION  ?= DOUBLE
PROFILE    ?= FALSE

DEBUG      ?= FALSE

DIM        ?= 2

COMP	   ?= gnu

USE_MPI    ?= FALSE
USE_OMP    ?= FALSE

NUMADV := 1

# This sets the EOS directory in Castro/EOS
EOS_DIR     ?= gamma_law

# This sets the network directory in Castro/Networks
NETWORK_DIR ?= general_null
GENERAL_NET_INPUTS ?= $(CASTRO_HOME)/Microphysics/networks/$(NETWORK_DIR)/gammalaw.net

Bpack   += $(TEST_DIR)/Make.package
Blocs   += $(TEST_DIR)/.

include $(CASTRO_HOME)/Exec/Make.Castro
