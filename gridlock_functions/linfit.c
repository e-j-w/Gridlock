long double confIntVal(long double input, const fit_results * fr, const data * d, int upper)
{
	int i;	
	long double cVal;
	long double maxCVal=-1.*BIG_NUMBER;
	long double minCVal=BIG_NUMBER;
	
	for (i=0;i<fr->ciEEValues;i++)
		{
			cVal=fr->ciEEVal[0][i]*input + fr->ciEEVal[1][i];
			//printf("i: %i, cVal: %LE\n",i,cVal);
			//printf("ciEEVals[%i]: %LE %LE\n",i,fr->ciEEVal[0][i],fr->ciEEVal[1][i]);
			if(cVal>maxCVal)
				maxCVal=cVal;
			if(cVal<minCVal)
				minCVal=cVal;
		}
	if(upper==1)
		return maxCVal;
	else
		return minCVal;

}

//prints the results
void printLin(const data * d, const parameters * p, fit_results * fr)
{
  //simplified data printing depending on verbosity setting
  if(p->verbose==1)
    {
      //print x and y intercept
      printf("%LE %LE\n",fr->fitVert[0],fr->fitVert[1]);
      return;
    }
  else if(p->verbose==2)
    {
      //print coefficient values
      int i;
      for(i=0;i<2;i++)
        printf("%LE ",fr->a[i]);
      printf("\n");
      return;
    }
  
  printf("\nFIT RESULTS\n-----------\n");
  printf("Uncertainties reported at 1-sigma.\n");
  printf("Fit function: f(x,y) = a1*x + a2\n\n");
  printf("Best chisq (fit): %0.3Lf\nBest chisq/NDF (fit): %0.3Lf\n\n",fr->chisq,fr->chisq/fr->ndf);
  printf("Coefficients from fit: a1 = %LE +/- %LE\n",fr->a[0],fr->aerr[0]);
  printf("                       a2 = %LE +/- %LE\n",fr->a[1],fr->aerr[1]);
  printf("\n");
  
  printf("x-intercept = %LE\n",fr->fitVert[0]);
  if ((float)(fr->fitVert[1]-confIntVal(0.0,fr,d,0))==(float)(confIntVal(0.0,fr,d,1)-fr->fitVert[1]))
    printf("y-intercept = %LE +/- %LE (from confidence interval)\n",fr->fitVert[1],fr->fitVert[1]-confIntVal(0.0,fr,d,0));
  else
    printf("y-intercept = %LE + %LE - %LE  (from confidence interval)\n",fr->fitVert[1],fr->fitVert[1]-confIntVal(0.0,fr,d,0),confIntVal(0.0,fr,d,1)-fr->fitVert[1]);
  
  /*//draw confidence interval ellipse
  handle=gnuplot_init();
  gnuplot_plot_xy(handle, fr->ciEEVal[0], fr->ciEEVal[1], fr->ciEEValues, "Data");
  getc(stdin);*/
    
}

void plotFormLin(const parameters * p, fit_results * fr)
{
	//set up equation forms for plotting
	if(strcmp(p->plotMode,"1d")==0)
		sprintf(fr->fitForm[0], "%LE*x + %LE",fr->a[0],fr->a[1]);
}

//fit data to a line of the form
//f(x) = a1*x + a2
void fitLin(const parameters * p, const data * d, fit_results * fr, plot_data * pd, int print)
{
  //construct equations
  int i,j;
  lin_eq_type linEq;
  linEq.dim=2;
  
  linEq.matrix[0][0]=d->xpowsum[0][2];
  linEq.matrix[0][1]=d->xpowsum[0][1]; 
  linEq.matrix[1][0]=d->xpowsum[0][1]; 
  linEq.matrix[1][1]=d->xpowsum[0][0];//bottom right entry
  
  linEq.vector[0]=d->mxpowsum[0][1];
  linEq.vector[1]=d->mxpowsum[0][0];
    
	//solve system of equations and assign values
	if(!(solve_lin_eq(&linEq)==1))
		{
			printf("ERROR: Could not determine fit parameters.\n");
			printf("Perhaps there are not enough data points to perform a fit?\n");
			exit(-1);
		}
  
  //save fit parameters  
  for(i=0;i<linEq.dim;i++)
    fr->a[i]=linEq.solution[i];
  long double f;
  fr->chisq=0;
  fr->ndf=d->lines-3;
  for(i=0;i<d->lines;i++)//loop over data points for chisq
    {
      f=fr->a[0]*d->x[0][i] + fr->a[1];
      fr->chisq+=(d->x[1][i] - f)*(d->x[1][i] - f)
                  /(d->x[1+1][i]*d->x[1+1][i]);
    }
  
  //Calculate covariances and uncertainties, see J. Wolberg 
  //'Data Analysis Using the Method of Least Squares' sec 2.5
  for(i=0;i<linEq.dim;i++)
    for(j=0;j<linEq.dim;j++)
      fr->covar[i][j]=linEq.inv_matrix[i][j]*(fr->chisq/fr->ndf);
  for(i=0;i<linEq.dim;i++)
    fr->aerr[i]=(long double)sqrt((double)(fr->covar[i][i]));
    
  //now that the fit is performed, use the fit parameters to find the intercept(s)
  fr->fitVert[0]=-1.0*fr->a[1]/fr->a[0];//x-intercept
  fr->fitVert[1]=fr->a[1];//y-intercept
  
  
	//generate slope/intercept pairs for confidence interval
	//Ref: A. Chester master thesis
  int gridSize=(int)(CI_EE_DIM/2.);
	long double c=2.30*fr->chisq/fr->ndf;//confidence level for 1-sigma in 2 parameters
	long double da0=sqrt((d->xpowsum[0][0]*c)/(d->xpowsum[0][0]*d->xpowsum[0][2] - d->xpowsum[0][1]*d->xpowsum[0][1]));
	//printf("sum: %LF\n",d->xpowsum[0][0]*d->xpowsum[0][2] - d->xpowsum[0][1]*d->xpowsum[0][1]);
	long double p0,p1,p2,a0b,a1b1,a1b2;
  p0=d->xpowsum[0][0];
  fr->ciEEValues=0;
	for (i=0;i<gridSize;i++)
		{
			a0b=fr->a[0]+((2.0*i - gridSize)/(double)gridSize)*da0;
			//printf("a0b: %Lf\n",a0b);
			p1=(2.*d->xpowsum[0][1]*(a0b - fr->a[0])) - 2.*d->xpowsum[0][0]*fr->a[1];
			p2=d->xpowsum[0][0]*fr->a[1]*fr->a[1] - 2.*d->xpowsum[0][1]*fr->a[1]*(a0b - fr->a[0]) + d->xpowsum[0][2]*(a0b - fr->a[0])*(a0b - fr->a[0]) - c;
			//printf("p1*p1 - 4.*p0*p2: %Lf\n",p1*p1 - 4.*p0*p2);
			a1b1=(-1.*p1 + sqrt(p1*p1 - 4.*p0*p2))/(2.*p0);
			a1b2=(-1.*p1 - sqrt(p1*p1 - 4.*p0*p2))/(2.*p0);
			//printf("i: %i, a0b: %LF, a1b1: %Lf, a1b2: %Lf\n",i,a0b,a1b1,a1b2);
			//record pairs of slope/intercept values on error ellipse
			fr->ciEEVal[0][fr->ciEEValues]=a0b;
			fr->ciEEVal[1][fr->ciEEValues]=a1b1;
			fr->ciEEValues++;
			fr->ciEEVal[0][fr->ciEEValues]=a0b;
			fr->ciEEVal[1][fr->ciEEValues]=a1b2;
			fr->ciEEValues++;
		}
	
	//construct the confidence interval
	for(i=0;i<CI_DIM;i++)
		{
			fr->ciXVal[0][i]=d->min_x[0] - (d->max_x[0]-d->min_x[0]) + (d->max_x[0]-d->min_x[0])*((3.0*i)/(CI_DIM-1.0));
			fr->ciUVal[0][i]=confIntVal(fr->ciXVal[0][i],fr,d,1);
			fr->ciLVal[0][i]=confIntVal(fr->ciXVal[0][i],fr,d,0);
			//printf("point %i, upper val: %lf, lower val: %lf\n",i,fr->ciUVal[0][i],fr->ciLVal[0][i]);
		}
	
	//print results
  if(print==1)
		printLin(d,p,fr);
	
	if((p->plotData==1)&&(p->verbose<1))
		{
			preparePlotData(d,p,fr,pd);
			plotFormLin(p,fr);
			plotData(p,fr,pd);
		}
  
}
