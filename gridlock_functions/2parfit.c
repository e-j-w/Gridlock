//fit data to a paraboloid of the form
//f(x,y) = a1*x^2 + a2*y^2 + a3*x*y + a4*x + a5*y + a6
void fit2Par(const data * d, fit_results * fr)
{
  //construct equations (n=2 specific case)
  int i,j;
  lin_eq_type linEq;
  linEq.dim=6;
  
  for(i=0;i<2;i++)//loop over free parameters
    for(j=i;j<2;j++)//loop over free parameters
      linEq.matrix[i][j]=d->xxpowsum[i][2][j][2];//top-left 2x2 entries
      
  linEq.matrix[0][2]=d->xxpowsum[0][3][1][1];
  linEq.matrix[0][3]=d->xpowsum[0][3];
  linEq.matrix[0][4]=d->xxpowsum[0][2][1][1];
  linEq.matrix[0][5]=d->xpowsum[0][2];
  
  linEq.matrix[1][2]=d->xxpowsum[0][1][1][3];
  linEq.matrix[1][3]=d->xxpowsum[0][1][1][2];
  linEq.matrix[1][4]=d->xpowsum[1][3];
  linEq.matrix[1][5]=d->xpowsum[1][2];
  
  linEq.matrix[2][2]=d->xxpowsum[0][2][1][2];
  linEq.matrix[2][3]=d->xxpowsum[0][2][1][1];
  linEq.matrix[2][4]=d->xxpowsum[0][1][1][2];
  linEq.matrix[2][5]=d->xxpowsum[0][1][1][1];
  
  linEq.matrix[3][3]=d->xpowsum[0][2];
  linEq.matrix[3][4]=d->xxpowsum[0][1][1][1];
  linEq.matrix[3][5]=d->xpowsum[0][1];
  
  linEq.matrix[4][4]=d->xpowsum[1][2];
  linEq.matrix[4][5]=d->xpowsum[1][1];
      
  linEq.matrix[5][5]=d->xpowsum[0][0];//bottom right entry
  
  //mirror the matrix (top right half mirrored to bottom left half)
  for(i=1;i<linEq.dim;i++)
    for(j=0;j<i;j++)
      linEq.matrix[i][j]=linEq.matrix[j][i];
  
  for(i=0;i<2;i++)
    linEq.vector[i]=d->mxpowsum[i][2];
  linEq.vector[2]=d->mxxpowsum[0][1][1][1];
  for(i=3;i<5;i++)
    linEq.vector[i]=d->mxpowsum[i-3][1];
  linEq.vector[5]=d->msum;
    
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
  fr->ndf=d->lines-7;
  for(i=0;i<d->lines;i++)//loop over data points to get chisq
    {
      f=fr->a[0]*d->x[0][i]*d->x[0][i] + fr->a[1]*d->x[1][i]*d->x[1][i] + fr->a[2]*d->x[0][i]*d->x[1][i] + fr->a[3]*d->x[0][i] + fr->a[4]*d->x[1][i] + fr->a[5];
      fr->chisq+=(d->x[2][i] - f)*(d->x[2][i] - f)/(d->x[2+1][i]*d->x[2+1][i]);;
    }
  //Calculate covariances and uncertainties, see J. Wolberg 
  //'Data Analysis Using the Method of Least Squares' sec 2.5
  for(i=0;i<linEq.dim;i++)
    for(j=0;j<linEq.dim;j++)
      fr->covar[i][j]=linEq.inv_matrix[i][j]*(fr->chisq/fr->ndf);
  for(i=0;i<linEq.dim;i++)
    fr->aerr[i]=(long double)sqrt((double)(fr->covar[i][i]));
    
  //now that the fit is performed, use the fit parameters (and the derivative of the fitting function) to find the vertex
  fr->fitVert[0]=4*fr->a[0]*fr->a[1]*fr->a[3] - 2*fr->a[0]*fr->a[2]*fr->a[4];
  fr->fitVert[0]/=2*fr->a[0]*fr->a[2]*fr->a[2] - 8*fr->a[0]*fr->a[0]*fr->a[1];
  fr->fitVert[1]=-1.*fr->a[2]*fr->fitVert[0] - fr->a[4];
  fr->fitVert[1]/=2*fr->a[1];

  //find the value of the fit function at the vertex
  fr->vertVal=fr->a[0]*fr->fitVert[0]*fr->fitVert[0] + fr->a[1]*fr->fitVert[1]*fr->fitVert[1] + fr->a[2]*fr->fitVert[0]*fr->fitVert[1] + fr->a[3]*fr->fitVert[0] + fr->a[4]*fr->fitVert[1] + fr->a[5];

  
}


//determine uncertainty ellipse bounds for the vertex by intersection of fit function with plane defining values at min + delta
//derived by: 
//1) setting f(x,y)=delta+min
//2) deriving x values as a function of y and vice versa via the quadratic formula
//3) setting the expression under the sqrts obtained to 0 to define bounds for x,y
//4) solving for upper and lower x,y bounds using the quadratic formula (calculated below)
void fit2ParChisqConf(fit_results * fr)
{
  
  long double a,b,c;
  long double delta=2.30;//confidence level for 1-sigma in 2 parameters
  fr->vertBoundsFound=1;
  
  a=4.*fr->a[1]*fr->a[0] - fr->a[2]*fr->a[2];
  b=4.*fr->a[1]*fr->a[3] - 2.*fr->a[2]*fr->a[4];
  c=4.*fr->a[1]*(fr->a[5] - delta - fr->vertVal) - fr->a[4]*fr->a[4];
  if((b*b - 4*a*c)<0.) 
    c=4.*fr->a[1]*(fr->a[5] + delta - fr->vertVal) - fr->a[4]*fr->a[4];//try flipping delta
  if((b*b - 4*a*c)<0.)  
    fr->vertBoundsFound=0;
  else
    {
      fr->vertUBound[0]=(-1.*b + (long double)sqrt((double)(b*b - 4*a*c)))/(2*a);
      fr->vertLBound[0]=(-1.*b - (long double)sqrt((double)(b*b - 4*a*c)))/(2*a);
    }
  
  a=4.*fr->a[0]*fr->a[1] - fr->a[2]*fr->a[2];
  b=4.*fr->a[0]*fr->a[4] - 2.*fr->a[2]*fr->a[3];
  c=4.*fr->a[0]*(fr->a[5] - delta - fr->vertVal) - fr->a[3]*fr->a[3];
  
  if((b*b - 4*a*c)<0.) 
    c=4.*fr->a[0]*(fr->a[5] + delta - fr->vertVal) - fr->a[3]*fr->a[3];//try flipping delta
  if((b*b - 4*a*c)<0.)  
    fr->vertBoundsFound=0;
  else
    {
      fr->vertUBound[1]=(-1.*b + (long double)sqrt((double)(b*b - 4*a*c)))/(2*a);
      fr->vertLBound[1]=(-1.*b - (long double)sqrt((double)(b*b - 4*a*c)))/(2*a);
    }

  //swap bounds if needed
  int i;
  for(i=0;i<2;i++)
    if(fr->vertLBound[i]>fr->vertUBound[i])
      {
        a=fr->vertUBound[i];
        fr->vertUBound[i]=fr->vertLBound[i];
        fr->vertLBound[i]=a;
      }

}


//prints fit data
void print2Par(const data * d, const parameters * p, const fit_results * fr)
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
  
  printf("\nFIT RESULTS\n-----------\n");
  printf("Uncertainties reported at 1-sigma.\n");
  printf("Fit function: f(x,y) = a1*x^2 + a2*y^2 + a3*x*y\n                     + a4*x + a5*y + a6\n\n");
  printf("Best chisq (fit): %0.3Lf\nBest chisq/NDF (fit): %0.3Lf\n\n",fr->chisq,fr->chisq/fr->ndf);
  printf("Coefficients from fit: a1 = %LE +/- %LE\n",fr->a[0],fr->aerr[0]);
  for(i=1;i<6;i++)
    printf("                       a%i = %LE +/- %LE\n",i+1,fr->a[i],fr->aerr[i]);
  printf("\n");
  
  if(fr->a[0]>=0)
    printf("Minimum in x direction, ");
  else
    printf("Maximum in x direction, ");
  if(fr->vertBoundsFound==1)
    {
      //these values were calculated at long double precision, 
      //check if they are the same to within float precision
      if ((float)(fr->fitVert[0]-fr->vertLBound[0])==(float)(fr->vertUBound[0]-fr->fitVert[0]))
        printf("x0 = %LE +/- %LE\n",fr->fitVert[0],fr->vertUBound[0]-fr->fitVert[0]);
      else
        printf("x0 = %LE + %LE - %LE\n",fr->fitVert[0],fr->vertUBound[0]-fr->fitVert[0],fr->fitVert[0]-fr->vertLBound[0]);
    }
  else
    printf("x0 = %LE\n",fr->fitVert[0]);
  if(fr->a[1]>=0)
    printf("Minimum in y direction, ");
  else
    printf("Maximum in y direction, ");
  if(fr->vertBoundsFound==1)
    {
      if ((float)(fr->fitVert[1]-fr->vertLBound[1])==(float)(fr->vertUBound[1]-fr->fitVert[1]))
        printf("y0 = %LE +/- %LE\n",fr->fitVert[1],fr->vertUBound[1]-fr->fitVert[1]);
      else
        printf("y0 = %LE + %LE - %LE\n",fr->fitVert[1],fr->vertUBound[1]-fr->fitVert[1],fr->fitVert[1]-fr->vertLBound[1]);
    }
  else
    printf("y0 = %LE\n",fr->fitVert[1]);
  
  printf("\nf(x0,y0) = %LE\n",fr->vertVal);
  
}

//generates the functional form of the fit function for plotting,
//which varies depending on the plotting mode (parameters may be fixed)
char * plotForm2Par(const parameters * p, const fit_results * fr, plot_data * pd, int plotNum)
{
  char * str;
  str=(char*)calloc(256,sizeof(char));
  if(strcmp(p->plotMode,"1d")==0)
    {
      if(plotNum==0)
        sprintf(str, "%Lf*(x**2) + %Lf*(%Lf**2) + %Lf*x*%Lf + %Lf*x + %Lf*%Lf + %Lf",fr->a[0],fr->a[1],pd->fixedParVal[1],fr->a[2],pd->fixedParVal[1],fr->a[3],fr->a[4],pd->fixedParVal[1],fr->a[5]);
      else if(plotNum==1)
        sprintf(str, "%Lf*(x**2) + %Lf*(%Lf**2) + %Lf*x*%Lf + %Lf*x + %Lf*%Lf + %Lf",fr->a[1],fr->a[0],pd->fixedParVal[0],fr->a[2],pd->fixedParVal[0],fr->a[4],fr->a[3],pd->fixedParVal[0],fr->a[5]);
    }
  else if(strcmp(p->plotMode,"2d")==0)
    {
      sprintf(str, "%Lf*(x**2) + %Lf*(y**2) + %Lf*x*y + %Lf*x + %Lf*y + %Lf",fr->a[0],fr->a[1],fr->a[2],fr->a[3],fr->a[4],fr->a[5]);
    }
  else
    {
      printf("ERROR: Invalid plot mode (%s), cannot get functional form.\n",p->plotMode);
      exit(-1);
    }
    
  return str;

}
