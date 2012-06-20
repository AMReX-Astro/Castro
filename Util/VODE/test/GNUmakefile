fsources = tvode.f

include ../GPackage.mak

F90 = g77
F90 = ifort -check all
F90 = lf95 -g --pca --nsav --ap --chk aesu --chkglobal

tvode: $(fsources)
	$(F90) -o tvode $^

vpath %.f . ../


clean:
	$(RM) *.o tvode

