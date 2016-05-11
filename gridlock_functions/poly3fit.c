//fit data to a 3rd order polynomial of the form
//f(x,y) = a1*x^3 + a2*x^2 + a3*x + a4
void fitPoly3(const data * d, fit_results * fr)
{
  //construct equations (n=1 specific case)
  int i,j;
  lin_eq_type linEq;
  linEq.dim=4;
  
  linEq.matrix[0][0]=d->xxpowsum[0][3][0][3];
  linEq.matrix[0][1]=d->xxpowsum[0][3][0][2];
  linEq.matrix[0][2]=d->xxpowsum[0][3][0][1];
  linEq.matrix[0][3]=d->xpowsum[0][3];
  
  linEq.matrix[1][1]=d->xxpowsum[0][3][0][1];
  linEq.matrix[1][2]=d->xpowsum[0][3];
  linEq.matrix[1][3]=d->xpowsum[0][2];
  
  linEq.matrix[2][2]=d->xpowsum[0][2];
  linEq.matrix[2][3]=d->xpowsum[0][1];
  
  linEq.matrix[3][3]=d->xpowsum[0][0];//bottom right entry
  
  //mirror the matrix (top right half mirrored to bottom left half)
  for(i=1;i<linEq.dim;i++)
    for(j=0;j<i;j++)
      linEq.matrix[i][j]=linEq.matrix[j][i];
  
  linEq.vector[0]=d->mxpowsum[0][3];
  linEq.vector[1]=d->mxpowsum[0][2];
  linEq.vector[2]=d->mxpowsum[0][1];
  linEq.vector[3]=d->mxpowsum[0][0];
    
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
  fr->ndf=d->lines-5;
  for(i=0;i<d->lines;i++)//loop over data points for chisq
    {
      f=fr->a[0]*d->x[0][i]*d->x[0][i]*d->x[0][i] + fr->a[1]*d->x[0][i]*d->x[0][i] + fr->a[2]*d->x[0][i] + fr->a[3];
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
    
  //now that the fit is performed, use the fit parameters (and the derivative of the fitting function) to find the critical points
  fr->fitVert[0]=-1.0*fr->a[1] - sqrt(fr->a[1]*fr->a[1] - 3.*fr->a[0]*fr->a[2]);
  fr->fitVert[0]/=3.*fr->a[0];
  fr->fitVert[1]=-1.0*fr->a[1] + sqrt(fr->a[1]*fr->a[1] - 3.*fr->a[0]*fr->a[2]);
  fr->fitVert[1]/=3.*fr->a[0];

  
}

//prints fit data
void printPoly3(const data * d, const parameters * p, const fit_results * fr)
{

  int i;

  //simplified data printing depending on verbosity setting
  if(p->verbose==1)
    {
      //print critical points
      printf("%LE %LE\n",fr->fitVert[0],fr->fitVert[1]);
      return;
    }
  
  printf("\nFIT RESULTS\n-----------\n");
  printf("Uncertainties reported at 1-sigma.\n");
  printf("Fit function: f(x,y) = a1*x^3 + a2*x^2 + a3*x + a4\n\n");
  printf("Best chisq (fit): %0.3Lf\nBest chisq/NDF (fit): %0.3Lf\n\n",fr->chisq,fr->chisq/fr->ndf);
  printf("Coefficients from fit: a1 = %LE +/- %LE\n",fr->a[0],fr->aerr[0]);
  for(i=1;i<4;i++)
    printf("                       a%i = %LE +/- %LE\n",i+1,fr->a[i],fr->aerr[i]);
  printf("\n");
  
  //check for NaN
  if((fr->fitVert[0]==fr->fitVert[0])&&(fr->fitVert[1]==fr->fitVert[1]))
    printf("Critical points at x = [ %LE %LE ]\n",fr->fitVert[0],fr->fitVert[1]);
  else
    printf("Fit function is monotonic (no critical points).\n");
  
}

//generates the functional form of the fit function for plotting
char * plotFormPoly3(const parameters * p, const fit_results * fr, plot_data * pd, int plotNum)
{
  char * str;
  str=(char*)calloc(256,sizeof(char));
  if(strcmp(p->plotMode,"1d")==0)
    {
      sprintf(str, "%Lf*(x**3) + %Lf*(x**2) + %Lf*x + %Lf",fr->a[0],fr->a[1],fr->a[2],fr->a[3]);
    }
  else
    {
      printf("ERROR: Invalid plot mode (%s), cannot get functional form.\n",p->plotMode);
      exit(-1);
    }
    
  return str;

}
