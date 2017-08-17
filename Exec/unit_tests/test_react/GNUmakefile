PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

DIM = 3

COMP = gnu

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE

USE_ACC = TRUE

# programs to be compiled
ALL: testburn.ex table

EOS_dir := helmholtz

Network_dir := aprox13

INTEGRATOR_DIR := BS

f90EXE_sources += testburn.f90

BLOCS = .
EXTERN_SEARCH = .

CASTRO_HOME := ../../..

include $(CASTRO_HOME)/Exec/Make.Castro


#my_objs = $(filter-out $(objEXETempDir)/main.o, $(objForExecs))

testburn.ex: $(objForExecs)
	@echo Linking $@ ...                                                    
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)

