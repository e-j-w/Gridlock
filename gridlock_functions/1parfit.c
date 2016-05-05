//fit data to a paraboloid of the form
//f(x,y) = a1*x^2 + a2*x + a3
void fit1Par(const data * d, fit_results * fr)
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
      printf("ERROR: Could not determine fit parameters.\n");
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
  
}

//determine uncertainty bounds for the vertex by intersection of fit function with line defining values at min + delta
//derived by: 
//1) setting f(x,y)=delta+min
//2) solving for x bounds using the quadratic formula (calculated below)
void fit1ParChisqConf(fit_results * fr)
{
  
  long double a,b,c;
  long double delta=1.00;//confidence level for 1-sigma in 1 parameter
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