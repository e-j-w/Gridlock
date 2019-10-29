//evaluates the fit function at the specified point
long double evalPoly4(long double x, const fit_results * fr)
{
	return fr->a[0]*x*x*x*x + fr->a[1]*x*x*x
					+ fr->a[2]*x*x + fr->a[3]*x + fr->a[4];
}

//evaluates fit function x values for the specified y value by iterative search
//y: value to search for
//closeToVal: where to start the search (in x)
//dir: search direction, 1 for increasing x, 0 for decreasing x
//resultVal: location to store the search result
//returns 1 if the search successful, 0 if not (fit function diverges past distanceAbortThreshold)
int evalPoly4XSearch(const data * d, const fit_results * fr, const long double y, const long double closeToVal, const int dir, long double * resultVal)
{

  //define step size in terms of data range
  long double stepSize = (d->max_x[0] - d->min_x[0])/10000.0;
  long double distanceAbortThreshold = (d->max_m - d->min_m)*20.0;
  long double distanceSolvedThreshold = (d->max_m - d->min_m)/10000000.0;
  long double searchVal=closeToVal;
  long double distance=evalPoly4(searchVal,fr)-y;
  int distSign;

  //printf("stepSize: %Lf, distanceSolvedThreshold: %Lf, distanceAbortThreshold: %Lfi\n",stepSize,distanceSolvedThreshold,distanceAbortThreshold);

  while((fabsl(distance))>distanceSolvedThreshold)
    {
      distSign = signbit(distance)==0;
      if(dir==1)
        searchVal+=stepSize;
      else
        searchVal-=stepSize;
      distance=evalPoly4(searchVal,fr)-y;
      if(distance>distanceAbortThreshold)
        return 0;
      if((signbit(distance)==0)!=distSign)
        {
          //we have passed the point of interest!
          //rewind back and decrease the step size
          if(dir==1)
            searchVal-=stepSize;
          else
            searchVal+=stepSize;
          stepSize/=10.;
        }
      //printf("searchVal: %Lf, distance: %Lf, abs(distance): %Lf\n",searchVal, distance,fabsl(distance));
    }

  memcpy(resultVal,&searchVal,sizeof(long double));
  //printf("Returning search result: %Lf\n",searchVal);
  return 1;

}


//determine uncertainty bounds for the local minimum by intersection of fit function with line defining values at min + delta
//done using iterative search
//delta is the desired confidence level (1.00 for 1-sigma in 1 parameter)
void fitPoly4ChisqConf(const parameters * p, const data * d, fit_results * fr, long double pt, int ind)
{
  
  long double vertVal = evalPoly4(pt,fr);
  //printf("vertval: %LF\n",vertVal);
  long double delta=p->ciDelta;

  long double boundVal=0.;
  fr->vertBoundsFound[ind]=0;

  if(evalPoly4XSearch(d, fr, vertVal+delta, pt, 1, &boundVal)==1)
    {
      memcpy(&fr->vertUBound[ind],&boundVal,sizeof(long double));
      if(evalPoly4XSearch(d, fr, vertVal+delta, pt, 0, &boundVal)==1)
        {
          memcpy(&fr->vertLBound[ind],&boundVal,sizeof(long double));
          fr->vertBoundsFound[ind]=1;
        }
    }
  
  //printf("Bounds: %Lf %Lf\n",fr->vertLBound[ind],fr->vertUBound[ind]);

}


//prints fit data
void printPoly4(const data * d, const parameters * p, const fit_results * fr)
{

  int i;

  //simplified data printing depending on verbosity setting
  if(p->verbose==1)
    {
      //print critical points
      printf("%LE %LE\n",fr->fitVert[0],fr->fitVert[1]);
      return;
    }
  else if(p->verbose==2)
    {
      //print coefficient values
      for(i=0;i<5;i++)
        printf("%LE ",fr->a[i]);
      printf("\n");
      return;
    } 
  
  printf("\nFIT RESULTS\n-----------\n");
  printf("Fit parameter uncertainties reported at 1-sigma.\n");
  printf("Fit function: f(x,y) = a1*x^4 + a2*x^3 + a3*x^2 + a4*x + a5\n\n");
  printf("Best chisq (fit): %0.3Lf\nBest chisq/NDF (fit): %0.3Lf\n\n",fr->chisq,fr->chisq/fr->ndf);
  printf("Coefficients from fit: a1 = %LE +/- %LE\n",fr->a[0],fr->aerr[0]);
  for(i=1;i<5;i++)
    printf("                       a%i = %LE +/- %LE\n",i+1,fr->a[i],fr->aerr[i]);
  printf("\n");
  
  printf("y-intercept = %LE\n",evalPoly4(0.0,fr));
  printf("\n");

  for(i=0;i<fr->numFitVert;i++)
    {
      //check that the vertex is a minimum
      if((12.0*fr->a[0]*fr->fitVert[i]*fr->fitVert[i] + 6.0*fr->a[1]*fr->fitVert[i] + 2.0*fr->a[2])>=0)
        {
          if(strcmp(p->dataType,"chisq")==0)
            {
              if(fr->vertBoundsFound[i]==1)
                {
                  if((float)(fr->vertUBound[i]-fr->fitVert[i])==(float)(fr->fitVert[i]-fr->vertLBound[i]))
                  printf("Local minimum (with %s confidence interval): x = %LE +/- %LE\n",p->ciSigmaDesc, fr->fitVert[i],fr->vertUBound[i]-fr->fitVert[i]);
                else
                  printf("Local minimum (with %s confidence interval): x = %LE + %LE - %LE\n",p->ciSigmaDesc,fr->fitVert[i],fr->vertUBound[i]-fr->fitVert[i],fr->fitVert[i]-fr->vertLBound[i]);
                }
              else
                {
                  printf("Local minimum: x = %LE\n",fr->fitVert[i]);
                  printf("Specified confidence interval (%s) is unbound for this local minimum.\n",p->ciSigmaDesc);
                }
            }
          else
            {
              printf("Local minimum: x = %LE\n",fr->fitVert[i]);
            }
        }
      else
        {
          printf("Local maximum: x = %LE\n",fr->fitVert[i]);
        }
        
    }
    
  
  if((p->findMinGridPoint == 1)||(p->findMaxGridPoint == 1)){
    printf("\n");
		int i;
    if(p->findMinGridPoint == 1){
      long double currentVal;
      long double minVal = BIG_NUMBER;
      int minPt = -1;
      for(i=0;i<d->lines;i++){
        currentVal = evalPoly4(d->x[0][i],fr);
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
        currentVal = evalPoly4(d->x[0][i],fr);
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


void plotFormPoly4(const parameters * p, fit_results * fr)
{
	//set up equation forms for plotting
	if(strcmp(p->plotMode,"1d")==0)
		{
			sprintf(fr->fitForm[0], "%.10LE*x*x*x*x + %0.20LE*x*x*x + %.10LE*x*x + %.10LE*x + %.10LE",fr->a[0],fr->a[1],fr->a[2],fr->a[3],fr->a[4]);
		}
}


//fit data to a 3rd order polynomial of the form
//f(x,y) = a1*x^4 + a2*x^3 + a3*x^2 + a4*x + a5
void fitPoly4(const parameters * p, const data * d, fit_results * fr, plot_data * pd, int print)
{

  int numFitPar = 5;
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
  linEq.dim=5;
  
  linEq.matrix[0][0]=d->xxpowsum[0][4][0][4];
  linEq.matrix[0][1]=d->xxpowsum[0][4][0][3];
  linEq.matrix[0][2]=d->xxpowsum[0][4][0][2];
  linEq.matrix[0][3]=d->xxpowsum[0][4][0][1];
  linEq.matrix[0][4]=d->xpowsum[0][4];
  
  linEq.matrix[1][1]=d->xxpowsum[0][3][0][3];
  linEq.matrix[1][2]=d->xxpowsum[0][3][0][2];
  linEq.matrix[1][3]=d->xxpowsum[0][3][0][1];
  linEq.matrix[1][4]=d->xpowsum[0][3];
  
  linEq.matrix[2][2]=d->xpowsum[0][4];
  linEq.matrix[2][3]=d->xpowsum[0][3];
  linEq.matrix[2][4]=d->xpowsum[0][2];

  linEq.matrix[3][3]=d->xpowsum[0][2];
  linEq.matrix[3][4]=d->xpowsum[0][1];
  
  linEq.matrix[4][4]=d->xpowsum[0][0];//bottom right entry
  
  //mirror the matrix (top right half mirrored to bottom left half)
  for(i=1;i<linEq.dim;i++)
    for(j=0;j<i;j++)
      linEq.matrix[i][j]=linEq.matrix[j][i];
  
  linEq.vector[0]=d->mxpowsum[0][4];
  linEq.vector[1]=d->mxpowsum[0][3];
  linEq.vector[2]=d->mxpowsum[0][2];
  linEq.vector[3]=d->mxpowsum[0][1];
  linEq.vector[4]=d->mxpowsum[0][0];
    
	//solve system of equations and assign values
	if(!(solve_lin_eq(&linEq)==1))
		{
			printf("ERROR: Could not determine fit parameters (poly4).\n");
			printf("Perhaps there are not enough data points to perform a fit?\n");
      printf("Otherwise you can also try adjusting the fit range using the UPPER_LIMITS and LOWER_LIMITS options.\n");
			exit(-1);
		}
  
  //save fit parameters  
  for(i=0;i<linEq.dim;i++)
    fr->a[i]=linEq.solution[i];
  long double f;
  fr->chisq=0;
  for(i=0;i<d->lines;i++)//loop over data points for chisq
    {
      f=fr->a[0]*d->x[0][i]*d->x[0][i]*d->x[0][i]*d->x[0][i] + fr->a[1]*d->x[0][i]*d->x[0][i]*d->x[0][i] + fr->a[2]*d->x[0][i]*d->x[0][i] + fr->a[3]*d->x[0][i] + fr->a[4];
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

  //find minima/maxima of fit
  //derivative of a quartic is a cubic, will use cubic functions (poly3fit.c)
  fit_results *svarfr=(fit_results*)calloc(1,sizeof(fit_results));
  svarfr->a[0]=4.0*fr->a[0];
  svarfr->a[1]=3.0*fr->a[1];
  svarfr->a[2]=2.0*fr->a[2];
  svarfr->a[3]=fr->a[3];
  //get cubic roots
  long double *roots=(long double*)calloc(3,sizeof(long double));
  fr->numFitVert=getPoly3Roots(0.0,svarfr,roots);
  for(i=0;i<fr->numFitVert;i++)
    fr->fitVert[i]=roots[i];
  free(roots);

  if(strcmp(p->dataType,"chisq")==0)
    {
      for(i=0;i<fr->numFitVert;i++)
        {
          //check that the vertex is a minimum
          if((12.0*fr->a[0]*fr->fitVert[i]*fr->fitVert[i] + 6.0*fr->a[1]*fr->fitVert[i] + 2.0*fr->a[2])>=0)
            {
              fitPoly4ChisqConf(p,d,fr,fr->fitVert[i],i);
            }
        }
    }

  //print results
  if(print==1)
		printPoly4(d,p,fr);
	
	if((p->plotData==1)&&(p->verbose<1))
		{
			preparePlotData(d,p,fr,pd);
			plotFormPoly4(p,fr);
			plotData(p,fr,pd);
		}
  
}
