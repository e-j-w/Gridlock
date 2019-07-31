//evaluates the fit function at the specified point
long double eval1Par(long double x, const fit_results * fr)
{
	return fr->a[0]*x*x+ fr->a[1]*x
					+ fr->a[2];
}

//evaluates the fit function x values at the specified y value
long double eval1ParX(long double y, const fit_results * fr, int pos)
{
  if(pos==0) //negative root
    return (-1.*fr->a[1] - sqrt(fr->a[1]*fr->a[1] - 4.*fr->a[0]*(fr->a[2]-y)))/(2.*fr->a[0]);
  else //positive root
    return (-1.*fr->a[1] + sqrt(fr->a[1]*fr->a[1] - 4.*fr->a[0]*(fr->a[2]-y)))/(2.*fr->a[0]);
}

//determine uncertainty bounds for the vertex by intersection of fit function with line defining values at min + delta
//delta is the desired confidence level (1.00 for 1-sigma in 1 parameter)
//derived by: 
//1) setting f(x,y)=delta+min
//2) solving for x bounds using the quadratic formula (calculated below)
void fit1ParChisqConf(const parameters * p, fit_results * fr)
{
  
  long double a,b,c;
  long double delta=p->ciDelta;
  fr->vertBoundsFound=1;
  
  a=fr->a[0];
  b=fr->a[1];
  c=fr->a[2] - delta - fr->vertVal;
  if((b*b - 4*a*c)<0.) 
    c=fr->a[2] + delta - fr->vertVal;//try flipping delta
  if((b*b - 4*a*c)<0.)
    fr->vertBoundsFound=0;
  else
    {
      fr->vertUBound[0]=(-1.*b + (long double)sqrt((double)(b*b - 4*a*c)))/(2*a);
      fr->vertLBound[0]=(-1.*b - (long double)sqrt((double)(b*b - 4*a*c)))/(2*a);
    }
  
  //swap bounds if needed
  if(fr->vertLBound[0]>fr->vertUBound[0])
    {
      a=fr->vertUBound[0];
      fr->vertUBound[0]=fr->vertLBound[0];
      fr->vertLBound[0]=a;
    }

}

//prints fit data
void print1Par(const data * d, const parameters * p, const fit_results * fr)
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
  printf("Fit function: f(x,y) = a1*x^2 + a2*x + a3\n\n");
  printf("Best chisq (fit): %0.3Lf\nBest chisq/NDF (fit): %0.3Lf\n\n",fr->chisq,fr->chisq/fr->ndf);
  printf("Coefficients from fit: a1 = %LE +/- %LE\n",fr->a[0],fr->aerr[0]);
  for(i=1;i<3;i++)
    printf("                       a%i = %LE +/- %LE\n",i+1,fr->a[i],fr->aerr[i]);
  printf("\n");
  
  if(fr->a[0]>=0)
    printf("Minimum in x direction");
  else
    printf("Maximum in x direction");
  if(fr->vertBoundsFound==1)
    {
      printf(" (with %s confidence interval), ",p->ciSigmaDesc);
      //these values were calculated at long double precision, 
      //check if they are the same to within float precision
      if ((float)(fr->fitVert[0]-fr->vertLBound[0])==(float)(fr->vertUBound[0]-fr->fitVert[0]))
        printf("x0 = %LE +/- %LE\n",fr->fitVert[0],fr->vertUBound[0]-fr->fitVert[0]);
      else
        printf("x0 = %LE + %LE - %LE\n",fr->fitVert[0],fr->vertUBound[0]-fr->fitVert[0],fr->fitVert[0]-fr->vertLBound[0]);
    }
  else
    {
      printf(", x0 = %LE\n",fr->fitVert[0]);
    }
    
  
  printf("\nf(x0) = %LE\n",fr->vertVal);

  printf("\ny-intercept = %LE\n",eval1Par(0.0,fr));

  if(strcmp(p->dataType,"chisq")==0)
    if( ((fr->a[0] >= 0.)&&(fr->fitVert[0]<0.)) || ((fr->a[0] < 0.)&&(fr->fitVert[0]>0.)) || (p->forceZeroX) )
      {
        if( ((fr->a[0] >= 0.)&&(fr->fitVert[0]<0.)) || ((fr->a[0] < 0.)&&(fr->fitVert[0]>0.)) )
          printf("Upper bound (with %s confidence interval) assuming minimum at zero: x = %LE\n",p->ciSigmaDesc, eval1ParX(eval1Par(0.0,fr)+p->ciDelta,fr,1));
        else
          printf("Lower bound (with %s confidence interval) assuming minimum at zero: x = %LE\n",p->ciSigmaDesc, eval1ParX(eval1Par(0.0,fr)+p->ciDelta,fr,0));
      }
      
      
      
  
  if((p->findMinGridPoint == 1)||(p->findMaxGridPoint == 1)){
    printf("\n");
    if(p->findMinGridPoint == 1){
      long double currentVal;
      long double minVal = BIG_NUMBER;
      int minPt = -1;
      for(i=0;i<d->lines;i++){
        currentVal = eval1Par(d->x[0][i],fr);
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
        currentVal = eval1Par(d->x[0][i],fr);
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

void plotForm1Par(const parameters * p, fit_results * fr)
{
	//set up equation forms for plotting
	if(strcmp(p->plotMode,"1d")==0)
		sprintf(fr->fitForm[0], "%.10LE*(x**2) + %.10LE*x + %.10LE",fr->a[0],fr->a[1],fr->a[2]);
}



//fit data to a paraboloid of the form
//f(x,y) = a1*x^2 + a2*x + a3
void fit1Par(const parameters * p, const data * d, fit_results * fr, plot_data * pd, int print)
{
  //construct equations (n=1 specific case)
  int i,j;
  lin_eq_type linEq;
  linEq.dim=3;
  
  linEq.matrix[0][0]=d->xpowsum[0][4];
  linEq.matrix[0][1]=d->xpowsum[0][3];
  linEq.matrix[0][2]=d->xpowsum[0][2];
  
  linEq.matrix[1][1]=d->xpowsum[0][2];
  linEq.matrix[1][2]=d->xpowsum[0][1];
  
  linEq.matrix[2][2]=d->xpowsum[0][0];//bottom right entry
  
  //mirror the matrix (top right half mirrored to bottom left half)
  for(i=1;i<linEq.dim;i++)
    for(j=0;j<i;j++)
      linEq.matrix[i][j]=linEq.matrix[j][i];
  
  linEq.vector[0]=d->mxpowsum[0][2];
  linEq.vector[1]=d->mxpowsum[0][1];
  linEq.vector[2]=d->mxpowsum[0][0];
    
	//solve system of equations and assign values
	if(!(solve_lin_eq(&linEq)==1))
		{
			printf("ERROR: Could not determine fit parameters (par1).\n");
			printf("Perhaps there are not enough data points to perform a fit?\n");
      printf("Otherwise you can also try adjusting the fit range using the UPPER_LIMITS and LOWER_LIMITS options.\n");
			exit(-1);
		}
  
  //save fit parameters  
  for(i=0;i<linEq.dim;i++)
    fr->a[i]=linEq.solution[i];
  long double f;
  fr->chisq=0;
  fr->ndf=d->lines-4;
  for(i=0;i<d->lines;i++)//loop over data points for chisq
    {
      f=fr->a[0]*d->x[0][i]*d->x[0][i] + fr->a[1]*d->x[0][i] + fr->a[2];
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
  fr->vertVal=fr->a[0]*fr->fitVert[0]*fr->fitVert[0] + fr->a[1]*fr->fitVert[0] + fr->a[2];
  
  
  if(strcmp(p->dataType,"chisq")==0)
  	fit1ParChisqConf(p,fr);//generate confidence interval bounds for chisq data
  
  //print results
  if(print==1)
		print1Par(d,p,fr);
	
	if((p->plotData==1)&&(p->verbose<1))
		{
			preparePlotData(d,p,fr,pd);
			plotForm1Par(p,fr);
			plotData(p,fr,pd);
		}
  
}
