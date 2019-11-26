MYBIN = /Users/russell/Lamont/PROJ_NoMelt_Cij/MINEOS_ACFLN/FORTRAN/bin
FC = gfortran
FFLAGS= $(MYFFLAGS)
# FFLAGS = -g -C
LFLAGS= -L$(MYLIB)


#-------------------------------------------------------------------------------
#f77 $(FFLAGS) $(LFLAGS) -o $(MYBIN)/mineos_qcorrectphv mineos_qcorrectphv.f
mineos_qcorrectphv: mineos_qcorrectphv.o 
	$(FC) $(FFLAGS) -o $(MYBIN)/mineos_qcorrectphv mineos_qcorrectphv.o
