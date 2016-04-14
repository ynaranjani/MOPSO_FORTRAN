# choose your compiler here:
FC = gfortran
#PGI compiler: pgf90
#gfortran: gfortran
#intel: ifort
FLAG = -O3 -g
#-O2 -stand f03  -check all -traceback -warn all -fstack-protector -assume protect_parens -implicitnone
#PGI flag: -mp -fast
#gfortran: -fopenmp -fbacktrace -fdump-parse-tree 
#intel: -openmp
OBJECTS= ProblemBank.o MOPSO.o MOPSOtools.o hyper.o shrink.o  rouletteWheel.o mutation.o
TARGET = main
TARGET: $(OBJECTS)
	$(FC) -o $(TARGET) $(FLAG) $(OBJECTS)
ProblemBank.mod: ProblemBank.o ProblemBank.f90
	$(FC)  -c $(FLAG) ProblemBank.f90
ProblemBank.o: ProblemBank.f90
	$(FC)  -c $(FLAG) ProblemBank.f90
MOPSO.o: ProblemBank.mod MOPSO.f90
	$(FC)  -c $(FLAG) MOPSO.f90
MOPSOtools.o: MOPSOtools.f90
	$(FC)  -c $(FLAG) MOPSOtools.f90
hyper.o: hyper.f90
	$(FC)  -c $(FLAG) hyper.f90	
shrink.o: shrink.f90
	$(FC)  -c $(FLAG) shrink.f90	
rouletteWheel.o: rouletteWheel.f90
	$(FC)  -c $(FLAG) rouletteWheel.f90	
mutation.o: mutation.f90
	$(FC)  -c $(FLAG) mutation.f90	

clean:
	rm *.o *.mod
	$(RM) $(TARGET)
