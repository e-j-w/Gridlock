CFLAGS 	= -I./lin_eq_solver -I./gridfit_functions -o2

all: lib gridfit

gridfit: gridfit.c gridfit.h lin_eq_solver.o
	gcc $(CFLAGS) gridfit.c -Wall -o gridfit lin_eq_solver.o
	
lib:lin_eq_solver/lin_eq_solver.c lin_eq_solver/lin_eq_solver.h
	gcc -I./lin_eq_solver -o2 -c -o lin_eq_solver.o lin_eq_solver/lin_eq_solver.c

clean:
	rm -rf *~ gridfit *.o
