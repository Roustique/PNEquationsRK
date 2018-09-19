G=gfortran
main.out: Runge-Kutta.o Main.o
	$(G) $^ -o $@
Main.o: Main.f95 my_prec.mod
	$(G) -c $<
Runge-Kutta.o: Runge-Kutta.f95 my_prec.mod
	$(G) -c $<
my_prec.mod: my_prec.f95
	$(G) -c $<
exec: main.out
	./main.out < Data.dat > result.dat
plot: Graphic.R
	./Graphic.R
