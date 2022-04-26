//forward declarations
void generateSums(data *,const parameters *);
void refitFilterPoly2Root0(const parameters *, const data *, fit_results *, plot_data *, long double);

//evaluates the fit function at the specified point
long double evalPoly2Root0(long double x, const fit_results * fr)
{
	return fr->a[0]*x*x+ fr->a[1]*x;
}

//evaluates the fit function x values at the specified y value
long double evalPoly2Root0X(long double y, const fit_results * fr, int pos)
{
  if(pos==0) //negative root
    return (-1.*fr->a[1] - sqrt(fr->a[1]*fr->a[1] - 4.*fr->a[0]*(0.-y)))/(2.*fr->a[0]);
  else //positive root
    return (-1.*fr->a[1] + sqrt(fr->a[1]*fr->a[1] - 4.*fr->a[0]*(0.-y)))/(2.*fr->a[0]);
}

//determine uncertainty bounds for the vertex by intersection of fit function with line defining values at min + delta
//delta is the desired confidence level (1.00 for 1-sigma in 1 parameter)
//derived by: 
//1) setting f(x,y)=delta+min
//2) solving for x bounds using the quadratic formula (calculated below)
void fitPoly2Root0ChisqConf(const parameters * p, fit_results * fr)
{
  
  long double a,b,c;
  long double delta=p->ciDelta;
  fr->vertBoundsFound[0]=1;
  
  a=fr->a[0];
  b=fr->a[1];
  c=0. - delta - fr->vertVal;
  if((b*b - 4*a*c)<0.) 
    c=delta - fr->vertVal;//try flipping delta
  if((b*b - 4*a*c)<0.)
    fr->vertBoundsFound[0]=0;
  else{
    fr->vertUBound[0]=(-1.*b + (long double)sqrt((double)(b*b - 4*a*c)))/(2*a);
    fr->vertLBound[0]=(-1.*b - (long double)sqrt((double)(b*b - 4*a*c)))/(2*a);
  }
  
  //swap bounds if needed
  if(fr->vertLBound[0]>fr->vertUBound[0]){
    a=fr->vertUBound[0];
    fr->vertUBound[0]=fr->vertLBound[0];
    fr->vertLBound[0]=a;
  }

}

//prints fit data
void printPoly2Root0(const data * d, const parameters * p, const fit_results * fr)
{

  int i;

  //simplified data printing depending on verbosity setting
  if(p->verbose==1)
    {
      //print vertex of paraboloid
      for(i=0;i<p->numVar;i++)
        printf("%LE ",fr->fitVert[i]);
      printf("\n");
      return;
    }
  else if(p->verbose==2)
    {
      //print coefficient values
      for(i=0;i<3;i++)
        printf("%LE ",fr->a[i]);
      printf("\n");
      return;
    } 
    
  printf("\nFIT RESULTS\n-----------\n");
  printf("Fit parameter uncertainties reported at 1-sigma.\n");
  printf("Fit function: f(x,y) = a1*x^2 + a2*x\n\n");
  printf("Best chisq (fit): %0.3Lf\nBest chisq/NDF (fit): %0.3Lf\n\n",fr->chisq,fr->chisq/fr->ndf);
  printf("Coefficients from fit: a1 = %LE +/- %LE\n",fr->a[0],fr->aerr[0]);
  printf("                       a2 = %LE +/- %LE\n",fr->a[1],fr->aerr[1]);
  printf("\n");
  
  if(fr->a[0]>=0)
    printf("Minimum in x direction");
  else
    printf("Maximum in x direction");
  if(fr->vertBoundsFound[0]==1){
    printf(" (with %s confidence interval): ",p->ciSigmaDesc);
    //these values were calculated at long double precision, 
    //check if they are the same to within float precision
    if ((float)(fr->fitVert[0]-fr->vertLBound[0])==(float)(fr->vertUBound[0]-fr->fitVert[0]))
      printf("x0 = %LE +/- %LE\n",fr->fitVert[0],fr->vertUBound[0]-fr->fitVert[0]);
    else
      printf("x0 = %LE + %LE - %LE\n",fr->fitVert[0],fr->vertUBound[0]-fr->fitVert[0],fr->fitVert[0]-fr->vertLBound[0]);
  }else{
    printf(": x0 = %LE\n",fr->fitVert[0]);
  }
  printf("f(x0) = %LE\n",fr->vertVal);

  printf("\ny-intercept is zero for this fit function.\n");

  if(strcmp(p->dataType,"chisq")==0)
    if( ((fr->a[0] >= 0.)&&(fr->fitVert[0]<0.)) || ((fr->a[0] < 0.)&&(fr->fitVert[0]>0.)) || (p->forceZeroX) )
      {
        if( ((fr->a[0] >= 0.)&&(fr->fitVert[0]<0.)) || ((fr->a[0] < 0.)&&(fr->fitVert[0]>0.)) || ((fr->a[0] >= 0.)&&(p->forceZeroX)) )
          {
            long double uBound = evalPoly2Root0X(evalPoly2Root0(0.0,fr)+p->ciDelta,fr,1);
            if(uBound == uBound)
              printf("Upper bound (with %s confidence interval) assuming minimum at zero: x = %LE\n",p->ciSigmaDesc, uBound);
            else
              printf("%s confidence interval is unbound if assuming a minimum at zero.\n", p->ciSigmaDesc);
          }
        else
          {
            long double lBound = evalPoly2Root0X(evalPoly2Root0(0.0,fr)-p->ciDelta,fr,0);
            if(lBound == lBound)
              printf("Lower bound (with %s confidence interval) assuming maximum at zero: x = %LE\n",p->ciSigmaDesc, lBound);
            else
              printf("%s confidence interval is unbound if assuming a maximum at zero.\n", p->ciSigmaDesc);
          }
          
      }
      
      
      
  
  if((p->findMinGridPoint == 1)||(p->findMaxGridPoint == 1)){
    printf("\n");
    if(p->findMinGridPoint == 1){
      long double currentVal;
      long double minVal = BIG_NUMBER;
      int minPt = -1;
      for(i=0;i<d->lines;i++){
        currentVal = evalPoly2Root0(d->x[0][i],fr);
        if(currentVal < minVal){
          minVal = currentVal;
          minPt = i;
        }
      }
      if(minPt >= 0){
        printf("Grid point corresponding to the lowest value (%LE) of the fitted function is at [ %0.3LE ].\n",minVal,d->x[0][minPt]);
      }
    }
    if(p->findMaxGridPoint == 1){
      long double currentVal;
      long double maxVal = -1.0*BIG_NUMBER;
      int maxPt = -1;
      for(i=0;i<d->lines;i++){
        currentVal = evalPoly2Root0(d->x[0][i],fr);
        if(currentVal > maxVal){
          maxVal = currentVal;
          maxPt = i;
        }
      }
      if(maxPt >= 0){
        printf("Grid point corresponding to the highest value (%LE) of the fitted function is at [ %0.3LE ].\n",maxVal,d->x[0][maxPt]);
      }
    }
  }

}

void plotFormPoly2Root0(const parameters * p, fit_results * fr)
{
	//set up equation forms for plotting
	if(strcmp(p->plotMode,"1d")==0)
		sprintf(fr->fitForm[0], "%.10LE*(x**2) + %.10LE*x",fr->a[0],fr->a[1]);
}



//fit data to a paraboloid of the form
//f(x,y) = a1*x^2 + a2*x
void fitPoly2Root0(const parameters * p, const data * d, fit_results * fr, plot_data * pd, int print)
{

  int numFitPar = 2;
  fr->ndf=d->lines-numFitPar;
  if(fr->ndf < 0)
    {
      printf("\nERROR: not enough data points for a fit (NDF < 0) using the %s function.\n",p->fitType);
      printf("%i data point(s) provided, %i data points needed.\n",d->lines,numFitPar);
      exit(-1);
    }
  else if(fr->ndf == 0)
    {
      if(p->verbose<1)
        {
          printf("\nWARNING: number of data points is equal to the number of fit parameters (%i).\n",numFitPar);
          printf("Fit is constrained to pass through data points (NDF = 0).\n");
        }
    }

  //construct equations (n=1 specific case)
  int i,j;
  lin_eq_type linEq;
  linEq.dim=2;
  
  linEq.matrix[0][0]=d->xpowsum[0][4];
  linEq.matrix[0][1]=d->xpowsum[0][3];
  
  linEq.matrix[1][1]=d->xpowsum[0][2];//bottom right entry
  linEq.matrix[1][0]=linEq.matrix[0][1];
  
  linEq.vector[0]=d->mxpowsum[0][2];
  linEq.vector[1]=d->mxpowsum[0][1];
    
	//solve system of equations and assign values
	if(!(solve_lin_eq(&linEq)==1))
		{
			printf("ERROR: Could not determine fit parameters (poly2root0).\n");
			printf("Perhaps there are not enough data points to perform a fit?\n");
      printf("Otherwise you can also try adjusting the fit range using the UPPER_LIMITS and LOWER_LIMITS options.\n");
			exit(-1);
		}
  
  //save fit parameters  
  for(i=0;i<linEq.dim;i++)
    fr->a[i]=linEq.solution[i];
  //refit filter  
  if(p->refitFilter==1){
    refitFilterPoly2Root0(p,d,fr,pd,p->refitFilterDist);
    return;
  }
  long double f;
  fr->chisq=0;
  for(i=0;i<d->lines;i++)//loop over data points for chisq
    {
      f=fr->a[0]*d->x[0][i]*d->x[0][i] + fr->a[1]*d->x[0][i];
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
    
  //now that the fit is performed, use the fit parameters (and the derivative of the fitting function) to find the minimum
  fr->fitVert[0]=-1.0*fr->a[1]/(2.*fr->a[0]);

  //find the value of the fit function at the vertex
  fr->vertVal=fr->a[0]*fr->fitVert[0]*fr->fitVert[0] + fr->a[1]*fr->fitVert[0];
  
  
  if(strcmp(p->dataType,"chisq")==0)
  	fitPoly2Root0ChisqConf(p,fr);//generate confidence interval bounds for chisq data
  
  //print results
  if(print==1)
		printPoly2Root0(d,p,fr);
	
	if((p->plotData==1)&&(p->verbose<1))
		{
			preparePlotData(d,p,fr,pd);
			plotFormPoly2Root0(p,fr);
			plotData(p,fr,pd);
		}
  
}

//generate a new data set and refit
void refitFilterPoly2Root0(const parameters * p, const data * d, fit_results * fr, plot_data * pd, long double distance){
  int i,j;

  //generate a new data set containing only filtered data
  data *nd=(data*)calloc(1,sizeof(data));
  parameters *np=(parameters*)calloc(1,sizeof(parameters));
  memcpy(np,p,sizeof(parameters));
  np->refitFilter=0;
  nd->lines=0;

  for (i=0;i<d->lines;i++){
    long double diff = fabs(d->x[p->numVar][i] - evalPoly2Root0(d->x[0][i],fr));
    if(diff<=distance){
      //keep this data point (remember to copy weight as well)
      for(j=0;j<=p->numVar+1;j++)
        nd->x[j][nd->lines] = d->x[j][i];
      nd->lines++;
    }
  }

  nd->max_m=d->max_m;
  nd->min_m=d->min_m;
  for(i=0;i<p->numVar;i++)
    nd->max_x[i]=d->max_x[i];

  if(p->verbose<1)
    printf("\nRefit filter: %i of %i data point(s) retained.\n",nd->lines,d->lines);
  //printDataInfo(nd,np); //see print_data_info.c

  generateSums(nd,np); //construct sums for fitting (see generate_sums.c)
  fitPoly2Root0(np,nd,fr,pd,1);

  free(nd);
  free(np);
}
