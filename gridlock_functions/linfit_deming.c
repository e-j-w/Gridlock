//fit data to a line of the form
//f(x) = a1*x + a2
//assuming errors in both x and y (Deming regression)
//see https://en.wikipedia.org/wiki/Deming_regression
//delta = ratio of variance in y/variance in x (1=errors perpendicular to line)
//reduces to regular linear fit for large delta
void fitLinDeming(const parameters * p, const data * d, fit_results * fr, long double delta)
{

  int i;
  long double xb,yb;
  long double sxx=0.;
  long double sxy=0.;
  long double syy=0.;
  
  xb=d->xpowsum[0][1]/d->lines;
  yb=d->mxpowsum[0][0]/d->lines;
  
  //compute 2nd degree sample moments
  for(i=0;i<d->lines;i++)
  	{
  		sxx+=(d->x[0][i] - xb)*(d->x[0][i] - xb);
  		sxy+=(d->x[0][i] - xb)*(d->x[1][i] - yb);
  		syy+=(d->x[1][i] - yb)*(d->x[1][i] - yb);
  	}
  sxx/=(d->lines + 1.0);
  sxy/=(d->lines + 1.0);
  syy/=(d->lines + 1.0);
  
  //compute fit parameters
  fr->a[0]=(syy - delta*sxx + sqrt((syy-delta*sxx)*(syy-delta*sxx) + 4.*delta*sxy*sxy))/(2.*sxy);
  fr->a[1]=yb - fr->a[0]*xb;
  
  
  long double x,y;
  fr->chisq=0;
  fr->ndf=d->lines-3;
  for(i=0;i<d->lines;i++)//loop over data points for chisq
    {
    	x=d->x[0][i] + (fr->a[0]/(fr->a[0]*fr->a[0] + delta))*(d->x[1][i] - fr->a[1] - fr->a[0]*d->x[0][i]);
    	y=fr->a[0]*x + fr->a[1];
      fr->chisq+=((d->x[1][i] - y)*(d->x[1][i] - y) + delta*(d->x[0][i] - x)*(d->x[0][i] - x))
                  /(d->x[1+1][i]*d->x[1+1][i]);
    }
  
  //Calculate covariances and uncertainties, see J. Wolberg 
  //'Data Analysis Using the Method of Least Squares' sec 2.5
  /*for(i=0;i<linEq.dim;i++)
    for(j=0;j<linEq.dim;j++)
      fr->covar[i][j]=linEq.inv_matrix[i][j]*(fr->chisq/fr->ndf);
  for(i=0;i<linEq.dim;i++)
    fr->aerr[i]=(long double)sqrt((double)(fr->covar[i][i]));*/
    
  //now that the fit is performed, use the fit parameters to find the intercept(s)
  fr->fitVert[0]=-1.0*fr->a[1]/fr->a[0];//x-intercept
  fr->fitVert[1]=fr->a[1];//y-intercept
  
  //set up equation forms for plotting
  if(strcmp(p->plotMode,"1d")==0)
    {
      sprintf(fr->fitForm[0], "%Lf*x + %Lf",fr->a[0],fr->a[1]);
    }
    
}

//prints the results
void printLinDeming(const data * d, const parameters * p, const fit_results * fr)
{
  //simplified data printing depending on verbosity setting
  if(p->verbose==1)
    {
      //print x and y intercept
      printf("%LE %LE\n",fr->fitVert[0],fr->fitVert[1]);
      return;
    }
  
  printf("\nFIT RESULTS\n-----------\n");
  //printf("Uncertainties reported at 1-sigma.\n");
  printf("Fit function: f(x,y) = a1*x + a2 (Deming regression, delta = %0.3Lf)\n\n",p->fitOpt);
  printf("Best chisq (fit): %0.3Lf\nBest chisq/NDF (fit): %0.3Lf\n\n",fr->chisq,fr->chisq/fr->ndf);
  printf("Coefficients from fit: a1 = %LE\n",fr->a[0]);
  printf("                       a2 = %LE\n",fr->a[1]);
  printf("\n");
  
  printf("x-intercept = %LE\n",fr->fitVert[0]);
  printf("y-intercept = %LE\n",fr->fitVert[1]);
  
  //printf("value at x=90 = %LE\n",fr->a[0]*90. + fr->a[1]);
  //printf("CI at x=90 = [%LE %LE]\n",confIntVal(90.,fr,d,1),confIntVal(90.,fr,d,0));
    
}
