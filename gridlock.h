#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "gnuplot_i.h"
#include "lin_eq_solver.h"

#define POWSIZE         12
#define MAXFILELENGTH   50000
#define BIG_NUMBER      1E10
#define NUM_LIST        5

typedef struct
{
  char filename[256];//name of the data file
  char fitType[256];//the type of fit (linear,parabola,etc)
  char dataType[256];//the type of data provided (regular, chisq values, etc.)
  int plotData;//0=don't plot, 1=plot
  char plotMode[256];//the plotting style to be used
  int numVar;
  long double ulimit[POWSIZE],llimit[POWSIZE];
  int verbose;//0=print everything,1=print vertex location only
  int readWeights;//0=don't read data weights from file,1=read weights
}parameters;

typedef struct
{
  int lines;//number of data points
  long double x[POWSIZE][MAXFILELENGTH];//array containing data points from the file, indexed by variable # then data point # 
  long double msum;//sum of measurements
  long double xpowsum[POWSIZE][POWSIZE];//sums of (x1)^0, (x1)^1, (x1)^2, etc. indexed first by variable # then by power #
  long double mxpowsum[POWSIZE][POWSIZE];//sums of m*(x1)^0, m*(x1)^1, m*(x1)^2, etc. indexed first by variable # then by power #
  long double mxxpowsum[POWSIZE][POWSIZE][POWSIZE][POWSIZE];//sums of m*(y1)*(x1)^0, m*(y1)*(x1)^1, m*(y1)*(x1)^2, etc. indexed by variable 1 #, variable 1 power #, variable 2 #, variable 2 power #
  long double xxpowsum[POWSIZE][POWSIZE][POWSIZE][POWSIZE];//sums of (y1), (y1)*(x1), (y1)*(x1)^2, etc. indexed by variable 1 #, variable 1 power #, variable 2 #, variable 2 power #
  long double xxxpowsum[POWSIZE][POWSIZE][POWSIZE][POWSIZE][POWSIZE][POWSIZE];//sums of (z1)*(y1)*(x1), (z1)*(y1)*(x1)^2, etc. indexed by variable 1 #, variable 1 power #, variable 2 #, variable 2 power #, variable 3 #, variable 3 power #
}data;

typedef struct
{
  long double fixedParVal[POWSIZE];//values to fix parameters at when plotting in less dimensions than the data provides
  double data[POWSIZE][POWSIZE][MAXFILELENGTH];//array containing data points to be plotted, indexed by plot # then variable # then data point #
  double fit[POWSIZE][POWSIZE][MAXFILELENGTH];//array containing fit data to be plotted, indexed by plot # then variable # then data point #
  int plotDataSize[POWSIZE];
  int numPlots;
  int axisLabelStyle[POWSIZE][POWSIZE];//0=normal,1=scientific notation
}plot_data;

typedef struct
{
  long double a[MAX_DIM]; //array holding parameters (desribing parboloid) from chisq minimization
  long double aerr[MAX_DIM]; //array holding uncertainties in parameters
  long double covar[MAX_DIM][MAX_DIM];//covariance between parameters (specified by the two indices)
  long double fitVert[POWSIZE]; //the vertex of the fit paraboloid
  long double vertUBound[POWSIZE],vertLBound[POWSIZE];//upper and lower bounds of the vertex
  long double vertVal; //value of the fit function at the vertex;
  long double chisq,ndf;
  int vertBoundsFound;
}fit_results;

//evil global variables
gnuplot_ctrl *handle;
int plotOpen;//1 if plots are being displayed, 0 otherwise
