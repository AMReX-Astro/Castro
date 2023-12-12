# Define the location of the CASTRO top directory,
# if not already defined by an environment variable.

CASTRO_HOME ?= ../../..

PRECISION   ?= DOUBLE
PROFILE     ?= FALSE

DEBUG       ?= FALSE

DIM         := 3

COMP	    ?= gnu

USE_MPI     ?= TRUE
USE_OMP     ?= FALSE

USE_GRAV     := TRUE
USE_ROTATION := TRUE

EOS_DIR     ?= helmholtz
NETWORK_DIR ?= aprox13

PROBLEM_DIR ?= ./

Bpack   += $(PROBLEM_DIR)/Make.package
Blocs   += $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
