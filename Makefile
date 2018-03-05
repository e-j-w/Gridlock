CFLAGS 	= -I./src/gnuplot_i -I./utils -I./src -o2

all: lib gridlock

gridlock: src/gridlock.c src/gridlock.h src/gnuplot_i.o src/lin_eq_solver.o
	@echo Making gridlock...
	gcc $(CFLAGS) src/gridlock.c -Wall -o gridlock src/gnuplot_i.o src/lin_eq_solver.o -lm
	@echo Tidying up...
	rm -rf *~ src/*.o
	
lib:src/gnuplot_i/gnuplot_i.c src/gnuplot_i/gnuplot_i.h src/lin_eq_solver.c src/lin_eq_solver.h
	@echo Making libraries...
	gcc -I./src/gnuplot_i -o2 -c -o src/gnuplot_i.o src/gnuplot_i/gnuplot_i.c
	gcc -I./src -o2 -c -o src/lin_eq_solver.o src/lin_eq_solver.c
	

clean:
	rm -rf *~ gridlock src/*.o
