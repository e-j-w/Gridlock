#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "gnuplot_i.h"
#include "lin_eq_solver.h"

#define POWSIZE         12
#define MAXFILELENGTH   500000
#define CI_EE_DIM				100 //# of data points to evaluate confidence interval error ellipse on
#define CI_DIM					100 //# of data points to use when plotting confidence interval
#define BIG_NUMBER      1E10
#define NUM_LIST        5
#define PI        			3.1415926535897932384626433832795028841971693993751

typedef struct
{
  char filename[256];//name of the data file
  char fitType[256];//the type of fit (linear,parabola,etc)
  long double fitOpt;//fit option value (eg. delta for Deming regression) 
  char dataType[256];//the type of data provided (regular, chisq values, etc.)
  long double ciDelta;//delta value for confidence interval calculation
  int plotData;//0=don't plot, 1=plot
  char plotMode[256];//the plotting style to be used
  int plotCI;//0=don't plot confidence interval, 1=plot it
  int numVar;
  long double ulimit[POWSIZE],llimit[POWSIZE];//upper and lower limits for variable values
  long double dllimit,dulimit;//upper and lower limits for data values
  int verbose;//0=print everything,1=print vertex location only
  int readWeights;//0=don't read data weights from file,1=read weights from file
  int uniWeight;//0=uniform weight not specified,1=uniform weight specified
  long double uniWeightVal;//value specified for uniform weights
  int filter;//0=no filter,1=linear filter
  double filterSigma;//sigma value to be used for filter
  int forceZeroX,forceZeroY,forceZeroZ;//whether or not to attempt forcing the fitted minimum to zero for x,y,z
  int numCIEvalPts; //number of points to evaluate the confidence interval bounds at (where applicable)
  long double CIEvalPts[100]; //array of x values at which to evaluate the confidence interval at
}parameters;

typedef struct
{
  int lines;//number of data points
  long double x[POWSIZE][MAXFILELENGTH];//array containing data points from the file, indexed by variable # then data point #
  long double max_x[POWSIZE],min_x[POWSIZE],max_m,min_m;//maximum and minimum values
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
  double max_m,min_m;//maximum and minimum values
  double ciData[POWSIZE][POWSIZE][MAXFILELENGTH];//array containing confidence interval data points to be plotted, indexed by plot # then variable # then data point #
  double fit[POWSIZE][POWSIZE][MAXFILELENGTH];//array containing fit data to be plotted, indexed by plot # then variable # then data point #
  int numFitPlotPts;//number of data points reserved for plotting fit data
  int numFitPtsPerVar;
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
  //confidence interval data
  int ciEEValues;//number of values on the error ellipse/ellipsoid used to get confidence interval
  double ciEEVal[POWSIZE][CI_EE_DIM];//values on the error ellipse/ellipsoid used to get confidence interval
  double ciUVal[POWSIZE][CI_DIM];//values on the upper confidence interval curve
  double ciLVal[POWSIZE][CI_DIM];//values on the lower confidence interval curve
  double ciXVal[POWSIZE][CI_DIM];//x values on the confidence interval curve
  //fit forms
  char fitForm[POWSIZE][256];//string containing form of the fitted equation
  /*char ciUForm[POWSIZE][256];//string containing form of the upper confidence interval
  char ciLForm[POWSIZE][256];//string containing form of the lower confidence interval
  char piUForm[POWSIZE][256];//string containing form of the upper prediction interval
  char piLForm[POWSIZE][256];//string containing form of the lower prediction interval*/
}fit_results;

//evil global variables
gnuplot_ctrl *handle;
int plotOpen;//1 if plots are being displayed, 0 otherwise
