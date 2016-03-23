#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "lin_eq_solver.h"

#define POWSIZE   12
#define MAXFILELENGTH   50000
#define BIG_NUMBER      1E100

//evil global variables
FILE *inp;
int i,j,k,l,m,n,o,p,status;
char str[256], eq[256];
double value[3];
long double x[POWSIZE][MAXFILELENGTH];//array containing data points from the file, indexed by variable # then data point #
long double msum;//sum of measurements
long double xpowsum[POWSIZE][POWSIZE];//sums of (x1)^0, (x1)^1, (x1)^2, etc. indexed first by variable # then by power #
long double mxpowsum[POWSIZE][POWSIZE];//sums of m*(x1)^0, m*(x1)^1, m*(x1)^2, etc. indexed first by variable # then by power #
long double mxxpowsum[POWSIZE][POWSIZE][POWSIZE][POWSIZE];//sums of m*(y1)*(x1)^0, m*(y1)*(x1)^1, m*(y1)*(x1)^2, etc. indexed by variable 1 #, variable 1 power #, variable 2 #, variable 2 power #
long double xxpowsum[POWSIZE][POWSIZE][POWSIZE][POWSIZE];//sums of (y1), (y1)*(x1), (y1)*(x1)^2, etc. indexed by variable 1 #, variable 1 power #, variable 2 #, variable 2 power #
long double xxxpowsum[POWSIZE][POWSIZE][POWSIZE][POWSIZE][POWSIZE][POWSIZE];//sums of (z1)*(y1)*(x1), (z1)*(y1)*(x1)^2, etc. indexed by variable 1 #, variable 1 power #, variable 2 #, variable 2 power #, variable 3 #, variable 3 power #

long double a[MAX_DIM]; //array holding parameters from chisq minimization

//linear equation solver
lin_eq_type linEq;
