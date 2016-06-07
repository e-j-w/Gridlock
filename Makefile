CFLAGS 	= -I./gnuplot_i -I./utils -I./gridlock_functions -o2

all: lib gridlock

gridlock: gridlock.c gridlock.h gnuplot_i.o lin_eq_solver.o tstat.o
	@echo Making gridlock...
	gcc $(CFLAGS) gridlock.c -Wall -o gridlock gnuplot_i.o tstat.o lin_eq_solver.o -lm
	@echo Tidying up...
	rm -rf *~ *.o
	
lib:gnuplot_i/gnuplot_i.c gnuplot_i/gnuplot_i.h utils/lin_eq_solver.c utils/lin_eq_solver.h
	@echo Making libraries...
	gcc -I./gnuplot_i -o2 -c -o gnuplot_i.o gnuplot_i/gnuplot_i.c
	gcc -I./utils -o2 -c -o lin_eq_solver.o utils/lin_eq_solver.c
	gcc -I./utils -o2 -c -o tstat.o utils/tstat.c
	

clean:
	rm -rf *~ gridlock *.o
