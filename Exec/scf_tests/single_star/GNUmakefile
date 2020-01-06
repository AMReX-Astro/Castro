# Define the location of the CASTRO top directory,
# if not already defined by an environment variable.

CASTRO_HOME := ../../..

# Location of this directory. Useful if
# you're trying to compile this from another location.

TEST_DIR = $(CASTRO_HOME)/Exec/scf_tests/single_star

PRECISION   ?= DOUBLE
PROFILE     ?= FALSE

DEBUG       ?= FALSE

DIM         := 3

COMP	    ?= gcc

USE_MPI     ?= FALSE
USE_OMP     ?= FALSE

USE_GRAV     := TRUE
USE_ROTATION := TRUE

EOS_DIR     ?= helmholtz
NETWORK_DIR ?= aprox13

USE_PROB_PARAMS = TRUE

Bpack   += $(TEST_DIR)/Make.package
Blocs   += $(TEST_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
