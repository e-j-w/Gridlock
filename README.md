# **Gridlock**

**A program to fit grids of data**

## Description

Gridlock is a program for fitting grids of data points with various functions and finding fit properties (eg. confidence intervals, intercepts, vertices).  Fitting routines are available for data with up to 3 free parameters.
 

## Features

* Fits data on a grid of up to 3 free variables.
* Uses fast non-iterative linear least-squares method.
* Automatically calculates and displays certain properties of fits (confidence intervals, intercepts, vertices).
* Various polynomial type fitting functions implemented.
* Simple ASCII data file format (compatible with gnuplot, excel, etc.).
* Plotting of fits is available via a built-in interface to gnuplot (the program will still compile and run if gnuplot is not present, but plotting will be unavailable).
* Lots of esoteric extra functionality added in an ad-hoc manner as the maintainer sees fit.
* Written in C because reasons.

## How to install

Use `make` to compile.  Optional data plotting requires `gnuplot` to be installed.

To run the program from anywhere, move the resulting `gridlock` executable to any directory under your `$PATH` environment variable.

Tested using gcc and GNU make on Ubuntu 14.04/16.04, Scientific Linux/CentOS 6, CentOS 7, and Arch Linux (as of June 2019).  The code may not build on compilers/platforms without support for 128-bit floating point (__float128) values, which are used to obtain higher precision in fits.  Otherwise the code is self-contained and should work on more or less any Linux distro.

## How to use

The program is run from the command line, with the only argument being the path to a file containing the data to be fitted:

```
gridlock /path/to/data/file
```

In addition to the data itself, data files must contain a line specifying the desired fit function, formatted 'FIT type' where 'type' is the name of the fit function (eg. 'par1', see *Available fitting functions* section below).  Example data files are included in the `sample` directory.

It is possible to specify grid fitting limits for each parameter in data files (see the *Options* section below and the `sample_3par.txt` file for an example).

It is also possible to plot data (using `gnuplot`) by including a line in the data file specifying a plotting style, see the *Options* section below.


## Available fitting functions

### Functions with 1 free parameter:

|**Name**|**Description**|**Form**|
|:---:|:---:|:---:|
|**lin** | line | f(x) = a<sub>1</sub>x + a<sub>2</sub>|
|**lin_deming** | line (Deming regression) | f(x) = a<sub>1</sub>x + a<sub>2</sub> with errors in x (see note below)|
|**poly2** | parabola (2nd order polynomial) | f(x) = a<sub>1</sub>x<sup>2</sup> + a<sub>2</sub>x + a<sub>3</sub>|
|**poly3** | cubic polynomial | f(x) = a<sub>1</sub>x<sup>3</sup> + a<sub>2</sub>x<sup>2</sup> + a<sub>3</sub>x + a<sub>4</sub>|
|**poly4** | quartic polynomial | f(x) = a<sub>1</sub>x<sup>4</sup> + a<sub>2</sub>x<sup>3</sup> + a<sub>3</sub>x<sup>2</sup> + a<sub>4</sub>x + a<sub>5</sub>|

NOTE: The **lin\_deming** function can take an optional parameter specifying the ratio of variances in y and x data.  Default value is 1 (errors perpendicular to line).  The parameter is specified on the fit function line (eg. 'FIT lin\_deming 3' for variance in y data 3x that of x data). 

### Functions with 2 free parameters:

|**Name**|**Description**|**Form**|
|:---:|:---:|:---:|
|**par2** | bivariate parabola (paraboloid) | f(x,y) = a<sub>1</sub>x<sup>2</sup> + a<sub>2</sub>y<sup>2</sup> + a<sub>3</sub>xy + a<sub>4</sub>x + a<sub>5</sub>y + a<sub>6</sub>|
|**2parpoly3** | bivariate cubic polynomial | f(x,y) = a<sub>1</sub>x<sup>3</sup> + a<sub>2</sub>y<sup>3</sup> + a<sub>3</sub>x<sup>2</sup>y + a<sub>4</sub>xy<sup>2</sup> + a<sub>5</sub>x<sup>2</sup> + a<sub>6</sub>y<sup>2</sup> +a<sub>7</sub>xy + a<sub>8</sub>x + a<sub>9</sub>y + a<sub>10</sub>|

### Functions with 3 free parameters:

|**Name**|**Description**|**Form**|
|:---:|:---:|:---:|
|**par3** | trivariate parabola | f(x,y,z) = a<sub>1</sub>x<sup>2</sup> + a<sub>2</sub>y<sup>2</sup> + a<sub>3</sub>z<sup>2</sup> + a<sub>4</sub>xy + a<sub>5</sub>xz + a<sub>6</sub>yz + a<sub>7</sub>x + a<sub>8</sub>y + a<sub>9</sub>z + a<sub>10</sub>|


## Options

Add these options (eg. 'PLOT 1d') as a single line anywhere in the data file to enable them:

|**Option**|**Effect**|
|:---:|:---:|
| PLOT 1d | Shows plot(s) in one variable, using fixed values for any other variables corresponding to the closest data points to the local minimum/maximum of the fit function.|
| PLOT 2d | Shows plot(s) in two variables (surface plot), using fixed values for any other variables corresponding to the closest data points to the local minimum/maximum of the fit function.|
| PLOT 3d | Shows plot in three variables (colour-coded heatmap plot) for the data, with the local minimum/maximum of the fit function marked on the map.|
| DATA_TYPE chisq | Tells the program that the data provided corresponds to chi-square goodness of fit statistic values computed on a grid for each of the free parameters.|
| PARAMETERS | If used, the program will only output the fit parameters (coordinates of the fit paraboloid vertex, or x and y intercept for a linear fit), which can be useful for interfacing the program with shell scripts.|
| COEFFICIENTS | If used, the program will only output the fit coefficients (a<sub>1</sub>, a<sub>2</sub>, ... ), which can be useful for interfacing the program with shell scripts.|
| FIND_MIN_GRID_POINT_FROM_FIT | If used, the program will report the grid point (out of the available data points) which corresponds to the lowest value of the fitted function.|
| FIND_MAX_GRID_POINT_FROM_FIT | If used, the program will report the grid point (out of the available data points) which corresponds to the highest value of the fitted function.|
| WEIGHTED | Use weights for the data points, specified in a column after the data values.|
| UNIFORM_WEIGHT value | Use a single weight value for all data points, without the need to put an extra column in the data file.|
| LINEAR_FILTER sigma | An outlier filtering option.  Before fitting, filter the data to emphasize prominent linear features, by only keeping data falling within sigma standard deviations of the mean x/y value.  Only applicable to data with one free parameter.|
| REFIT_FILTER value | An outlier filtering option.  After performing the initial fit, drop all data which is a distance greater than 'value' away from the corresponding fit value, and then refit the data.|
| IGNORE_PAR par | Ignore a certain parameter ('par' may be 'x', 'y', or 'z') in the data when fitting and plotting.  Equivalent to removing the column of values corresponding to the specified parameter from the data file.|
| SLICE_PAR par value | Slice the grid at the specified parameter and value ('par' may be 'x', 'y', or 'z') in the data when fitting and plotting (ie. take only the data where the specified parameter has the specified value, and fit only the remaining parameters).|
| LOWER_LIMITS value1 value2 value3 | Lower fit limits for each variable (specify as many values as there are variables).  Use with UPPER_LIMITS to specify a fit range.|
| UPPER_LIMITS value1 value2 value3 | Upper fit limits for each variable (specify as many values as there are variables).  Use with LOWER_LIMITS to specify a fit range.|
| DATA_LOWER_LIMIT value | Lower fit limit for data values.  Use with DATA_UPPER_LIMIT to specify a fit range for data values.  This option can be used to filter outlier data.|
| DATA_UPPER_LIMIT value | Upper fit limit for data values.  Use with DATA_LOWER_LIMIT to specify a fit range for data values.  This option can be used to filter outlier data.|
| ZEROX | When fitting using the *par1*, *par2*, *poly3*, or *2parpoly3* functions with chisq data, show the fit result assuming a minimum in x at x=0.  This option can be stacked with ZEROY for the *par2* and *2parpoly3* fit functions.|
| ZEROY | When fitting using the *par2* or *2parpoly3* functions with chisq data, show the fit result assuming a minimum in y at y=0.  This option can be stacked with ZEROX.|
| EVAL_CI value | When fitting a function (such as 'lin') which provides a confidence interval, evaluate the bounds of the confidence interval for the given value of the independent variable.|
| SET_CI_SIGMA value | When fitting chi-square data, manually set the sigma value used to evaluate uncertainties (valid values are 1, 2, 3, 90%).  The program will then handle the appropriate confidence bounds for the number of free parameters used.  Default value is 1-sigma.|
| SET_CI_DELTA value | For people who know what they're doing and for whom SET_CI_SIGMA isn't enough.  When fitting chi-square data, manually set the delta value used to evaluate confidence bounds (by default, delta is set to the 1-sigma bound ie. 1.00 for 1 parameter, 2.30 for 2 parameters, etc.).|


## Acknowledgments

This program uses the public domain `gnuplot_i` library by N. Devillard for displaying plots.
