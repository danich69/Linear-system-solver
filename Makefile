Gauss: prog
	./prog Gauss
ModGauss: prog
	./prog ModGauss
Jordan: prog
	./prog Jordan
prog: main.o Linear_system.o
	ifort $^ -o $@ -debug
main.o: main.f90 linear_system.mod
	ifort $^ -c
linear_system.mod Linear_system.o: Linear_system.f90
	ifort Linear_system.f90 -c
clean: 
	rm -f *.o *mod cr prog
