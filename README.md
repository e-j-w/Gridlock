# **Gridlock** - a program to fit grids of data

Maintainer: Jonathan Williams


## Description

Gridlock is a program for fitting grids of data points with various functions and finding fit properties (eg. confidence intervals, intercepts, vertices).  Fitting routines are available for data with up to 3 free parameters.  Example data is included in the `sample_*.txt` files.

Data files must contain a line specifying the fit function, formatted 'FIT type' where 'type' is the name of the fit function (eg. 'par1', see *Fitting Functions Available* section below).  It is possible to specify grid fitting limits for each parameter in data files (see the *Options* section below and the `sample_3par.txt` file for an example).  It is also possible to plot data (using `gnuplot`) by including a line in the data file specifying a plotting style, see the options section below.
 

## Features

* Fits data on a grid of up to 3 free variables
* Uses fast non-iterative linear least-squares method
* Automatically calculates and displays certain properties of fits (confidence intervals, intercepts, vertices)
* Various polynomial type fitting functions implemented.
* Simple ASCII data file format (compatible with gnuplot, excel, etc.)
* Plotting of fits is available via a built-in interface to gnuplot (program will still compile and run if gnuplot is not present, but plotting will be unavailable)


## How to Install

Use `make` to compile.  Optional data plotting requires `gnuplot` to be installed.

To run the program from anywhere, move the resulting `gridlock` executable to any directory under your `$PATH` environment variable.

Tested using gcc and GNU make on Ubuntu 14.04/16.04, Scientific Linux/CentOS 6, and Arch Linux (as of July 2016).  The code is self-contained and should work on more or less any Linux distro.


## Fitting Functions Available

### Functions with 1 free parameter:

|**Name**|**Description**|**Form**|
|:---:|:---:|:---:|
|**lin** | line | f(x) = a1\*x + a2|
|**lin_deming** | line (Deming regression) | f(x) = a1\*x + a2 with errors in x (see note below)|
|**par1** | parabola (2nd order polynomial) | f(x) = a1\*x^2 + a2\*x + a3|
|**poly3** | cubic polynomial | f(x) = a1\*x^3 + a2\*x^2 + a3\*x + a4|

NOTE: The **lin\_deming** function can take an optional parameter specifying the ratio of variances in y and x data.  Default value is 1 (errors perpendicular to line).  The parameter is specified on the fit fucntion line (eg. 'FIT lin\_deming 3' for variance in y data 3x that of x data). 

### Functions with 2 free parameters:

|**Name**|**Description**|**Form**|
|:---:|:---:|:---:|
|**par2** | bivariate parabola (paraboloid) | f(x,y) = a1\*x^2 + a2\*y^2 + a3\*x\*y + a4\*x + a5\*y + a6|
|**2parpoly3** | bivariate cubic polynomial | f(x,y) = a1\*x^3 + a2\*y^3 + a3\*x^2\*y + a4\*x\*y^2 + a5\*x^2 + a6\*y^2 +a7\*x\*y + a8\*x + a9\*y + a10|

### Functions with 3 free parameters:

|**Name**|**Description**|**Form**|
|:---:|:---:|:---:|
|**par3** | trivariate parabola | f(x,y,z) = a1\*x^2 + a2\*y^2 + a3\*z^2 + a4\*x\*y + a5\*x\*z + a6\*y\*z + a7\*x + a8\*y + a9\*z + a10|



## Options

Add these options (eg. 'PLOT 1d') as a single line anywhere in the data file to enable them:

|**Option**|**Effect**|
|:---:|:---:|
| DATA_TYPE chisq | Tells the program that the data provided corresponds to reduced chi-square goodness of fit statistic values computed on a grid for each of the free parameters.|
| PARAMETERS | If used, the program will only output the fit parameters (coordinates of the fit paraboloid vertex, or x and y intercept for a linear fit), which can be useful for interfacing the program with shell scripts.|
| COEFFICIENTS | If used, the program will only output the fit coefficients (a1, a2, ... ), which can be useful for interfacing the program with shell scripts.|
| PLOT 1d | Shows plot(s) in one variable, using fixed values for any other variables corresponding to the closest data points to the local minimum/maximum of the fit function.|
| PLOT 2d | Shows plot(s) in two variables (surface plot), using fixed values for any other variables corresponding to the closest data points to the local minimum/maximum of the fit function.|
| PLOT 3d | Shows plot in three variables (colour-coded heatmap plot) for the data, with the local minimum/maximum of the fit function marked on the map.|
| WEIGHTED | Use weights for the data points, specified in a column after the data values.|
| UNIFORM_WEIGHT value | Use a single weight value for all data points, without the need to put an extra column in the data file.|
| LINEAR_FILTER sigma | Before fitting, filter the data to emphasize prominent linear features, by only keeping data falling within sigma standard deviations of the mean x/y value.  Only applicable to data with one free parameter.|
| LOWER_LIMITS value1 value2 value3 | Lower fit limits for each variable (specify as many values as there are variables).  Use with UPPER_LIMITS to specify a fit range.|
| UPPER_LIMITS value1 value2 value3 | Upper fit limits for each variable (specify as many values as there are variables).  Use with LOWER_LIMITS to specify a fit range.|
| DATA_LOWER_LIMIT value | Lower fit limit for data values.  Use with DATA_UPPER_LIMIT to specify a fit range for data values.|
| DATA_UPPER_LIMIT value | Upper fit limit for data values.  Use with DATA_LOWER_LIMIT to specify a fit range for data values.|


## Acknowledgments

This program uses the public domain `gnuplot_i` library by N. Devillard for displaying plots.
