FBIN = /Users/russell/Lamont/PROJ_NoMelt_Cij/MINEOS_ACFLN/FORTRAN/bin
FC = gfortran
FFLAGS=-ffixed-line-length-none 
#-L/usr/local/include 
#FFLAGS2=-march=x86_64

all:  $(FBIN)/frechet_ACFLN_love 

.f.o: 
	$(FC) $(FFLAGS) $(FFLAGS2) -c $*.f

#----------------------------------

$(FBIN)/frechet_ACFLN_love: frechet_ACFLN_love.f
	$(FC) $(FFLAGS) -o $(FBIN)/frechet_ACFLN_love frechet_ACFLN_love.f
	
clean: 
	rm -rf *.o $(FBIN)/frechet_ACFLN_love
