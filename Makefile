FC=gfortran
FFLAGS= -fcheck=all -ffpe-trap=invalid,zero,overflow -O0 -fbacktrace -g
SRC=module.f90 maxwell-in-fortran.f90
OBJ=$(SRC:.f90=.o)

%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

run: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

ex: run
	./run

plot: ex
	python Visualization.py

clean:
	rm -f figures/*.png *.o *.mod *.txt run 
