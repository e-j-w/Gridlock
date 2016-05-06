//fit data to a line of the form
//f(x,y) = a1*x + a2
void fitLin(const data * d, fit_results * fr)
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
  
}
