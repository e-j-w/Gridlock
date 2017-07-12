//evaluates the fit function at the specified point
long double eval3Par(long double x,long double y,long double z, const fit_results * fr)
{
	return fr->a[0]*x*x + fr->a[1]*y*y + fr->a[2]*z*z 
	+ fr->a[3]*x*y + fr->a[4]*x*z + fr->a[5]*y*z 
	+ fr->a[6]*x + fr->a[7]*y + fr->a[8]*z 
	+ fr->a[9];
}


//determine uncertainty ellipsoid bounds for the vertex by from fit function values fixed at min + delta
//derived using the same procedure as for 2 free variables (see 2parfit.c), with an extra step solving the quadratic formula in between
//delta is the desired confidence level (3.53 for 1-sigma in 3 parameters)
void fit3ParChisqConf(const parameters * p, fit_results * fr)
{
  
  long double a,b,c;
  long double delta=p->ciDelta;
  fr->vertBoundsFound=1;
  
  a=(16.*fr->a[1]*fr->a[2] - 4.*fr->a[5]*fr->a[5])*(4.*fr->a[0]*fr->a[2] - fr->a[4]*fr->a[4]) - 16.*(fr->a[2]*fr->a[2]*fr->a[3]*fr->a[3] - fr->a[2]*fr->a[3]*fr->a[4]*fr->a[5]) - 4.*fr->a[4]*fr->a[4]*fr->a[5]*fr->a[5];
  b=(16.*fr->a[1]*fr->a[2] - 4.*fr->a[5]*fr->a[5])*(4.*fr->a[2]*fr->a[6] - 2.*fr->a[4]*fr->a[8]) - 16.*(2.*fr->a[2]*fr->a[2]*fr->a[3]*fr->a[7] - fr->a[2]*fr->a[5]*(fr->a[3]*fr->a[8] + fr->a[7]*fr->a[4])) - 8.*fr->a[5]*fr->a[5]*fr->a[4]*fr->a[8];
  c=(16.*fr->a[1]*fr->a[2] - 4.*fr->a[5]*fr->a[5])*(4.*fr->a[2]*(fr->a[9] - delta - fr->vertVal) - fr->a[8]*fr->a[8]) - 16.*(fr->a[2]*fr->a[2]*fr->a[7]*fr->a[7] - fr->a[2]*fr->a[5]*fr->a[7]*fr->a[8]) - 4.*fr->a[5]*fr->a[5]*fr->a[8]*fr->a[8]; 
  if((b*b - 4*a*c)<0.)
    c=(16.*fr->a[1]*fr->a[2] - 4.*fr->a[5]*fr->a[5])*(4.*fr->a[2]*(fr->a[9] + delta - fr->vertVal) - fr->a[8]*fr->a[8]) - 16.*(fr->a[2]*fr->a[2]*fr->a[7]*fr->a[7] - fr->a[2]*fr->a[5]*fr->a[7]*fr->a[8]) - 4.*fr->a[5]*fr->a[5]*fr->a[8]*fr->a[8];//try flipping delta
  if((b*b - 4*a*c)<0.)
    fr->vertBoundsFound=0;
  else
    {
      fr->vertUBound[0]=(-1.*b + (long double)sqrt((double)(b*b - 4.*a*c)))/(2.*a);
      fr->vertLBound[0]=(-1.*b - (long double)sqrt((double)(b*b - 4.*a*c)))/(2.*a);
    }

  a=(16.*fr->a[0]*fr->a[2] - 4.*fr->a[4]*fr->a[4])*(4.*fr->a[0]*fr->a[1] - fr->a[3]*fr->a[3]) - 16.*(fr->a[0]*fr->a[0]*fr->a[5]*fr->a[5] - fr->a[0]*fr->a[3]*fr->a[4]*fr->a[5]) - 4.*fr->a[3]*fr->a[3]*fr->a[4]*fr->a[4];  
  b=(16.*fr->a[0]*fr->a[2] - 4.*fr->a[4]*fr->a[4])*(4.*fr->a[0]*fr->a[7] - 2.*fr->a[3]*fr->a[6]) - 16.*(2.*fr->a[5]*fr->a[8]*fr->a[0]*fr->a[0] - fr->a[0]*fr->a[4]*(fr->a[5]*fr->a[6] + fr->a[8]*fr->a[3])) - 8.*fr->a[4]*fr->a[4]*fr->a[3]*fr->a[6];  
  c=(16.*fr->a[0]*fr->a[2] - 4.*fr->a[4]*fr->a[4])*(4.*fr->a[0]*(fr->a[9] - delta - fr->vertVal) - fr->a[6]*fr->a[6]) - 16.*(fr->a[0]*fr->a[0]*fr->a[8]*fr->a[8] - fr->a[0]*fr->a[4]*fr->a[6]*fr->a[8]) - 4.*fr->a[4]*fr->a[4]*fr->a[6]*fr->a[6];
  if((b*b - 4*a*c)<0.) 
    c=(16.*fr->a[0]*fr->a[2] - 4.*fr->a[4]*fr->a[4])*(4.*fr->a[0]*(fr->a[9] + delta - fr->vertVal) - fr->a[6]*fr->a[6]) - 16.*(fr->a[0]*fr->a[0]*fr->a[8]*fr->a[8] - fr->a[0]*fr->a[4]*fr->a[6]*fr->a[8]) - 4.*fr->a[4]*fr->a[4]*fr->a[6]*fr->a[6];//try flipping delta
  if((b*b - 4*a*c)<0.)  
    fr->vertBoundsFound=0;
  else
    {
      fr->vertUBound[1]=(-1.*b + (long double)sqrt((double)(b*b - 4.*a*c)))/(2.*a);
      fr->vertLBound[1]=(-1.*b - (long double)sqrt((double)(b*b - 4.*a*c)))/(2.*a);
    }
  
  a=(16.*fr->a[0]*fr->a[1] - 4.*fr->a[3]*fr->a[3])*(4.*fr->a[0]*fr->a[2] - fr->a[4]*fr->a[4]) - 16.*(fr->a[0]*fr->a[0]*fr->a[5]*fr->a[5] - fr->a[0]*fr->a[3]*fr->a[4]*fr->a[5]) - 4.*fr->a[3]*fr->a[3]*fr->a[4]*fr->a[4];
  b=(16.*fr->a[0]*fr->a[1] - 4.*fr->a[3]*fr->a[3])*(4.*fr->a[0]*fr->a[8] - 2.*fr->a[4]*fr->a[6]) - 16.*(2.*fr->a[5]*fr->a[7]*fr->a[0]*fr->a[0] - fr->a[0]*fr->a[3]*(fr->a[5]*fr->a[6] + fr->a[7]*fr->a[4])) - 8.*fr->a[3]*fr->a[3]*fr->a[4]*fr->a[6];
  c=(16.*fr->a[0]*fr->a[1] - 4.*fr->a[3]*fr->a[3])*(4.*fr->a[0]*(fr->a[9] - delta - fr->vertVal) - fr->a[6]*fr->a[6]) - 16.*(fr->a[0]*fr->a[0]*fr->a[7]*fr->a[7] - fr->a[0]*fr->a[3]*fr->a[6]*fr->a[7]) - 4.*fr->a[3]*fr->a[3]*fr->a[6]*fr->a[6];  
  if((b*b - 4*a*c)<0.) 
    c=(16.*fr->a[0]*fr->a[1] - 4.*fr->a[3]*fr->a[3])*(4.*fr->a[0]*(fr->a[9] + delta - fr->vertVal) - fr->a[6]*fr->a[6]) - 16.*(fr->a[0]*fr->a[0]*fr->a[7]*fr->a[7] - fr->a[0]*fr->a[3]*fr->a[6]*fr->a[7]) - 4.*fr->a[3]*fr->a[3]*fr->a[6]*fr->a[6];//try flipping delta
  if((b*b - 4*a*c)<0.)  
    fr->vertBoundsFound=0;
  else
    {
      fr->vertUBound[2]=(-1.*b + (long double)sqrt((double)(b*b - 4.*a*c)))/(2.*a);
      fr->vertLBound[2]=(-1.*b - (long double)sqrt((double)(b*b - 4.*a*c)))/(2.*a);
    }
  
  //swap bounds if needed
  int i;
  for(i=0;i<3;i++)
    if(fr->vertLBound[i]>fr->vertUBound[i])
      {
        a=fr->vertUBound[i];
        fr->vertUBound[i]=fr->vertLBound[i];
        fr->vertLBound[i]=a;
      }
    

}

//prints fit data
void print3Par(const data * d, const parameters * p, const fit_results * fr)
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
      for(i=0;i<10;i++)
        printf("%LE ",fr->a[i]);
      printf("\n");
      return;
    }
  
  printf("\nFIT RESULTS\n-----------\n");
  printf("Uncertainties reported at 1-sigma.\n");
  printf("Fit function: f(x,y,z) = a1*x^2 + a2*y^2 + a3*z^2\n                       + a4*x*y + a5*x*z + a6*y*z\n                       + a7*x + a8*y + a9*z + a10\n\n");
  printf("Best chisq (fit): %0.3Lf\nBest chisq/NDF (fit): %0.3Lf\n\n",fr->chisq,fr->chisq/fr->ndf);
  printf("Coefficients from fit: a1 = %LE +/- %LE\n",fr->a[0],fr->aerr[0]);
  for(i=1;i<10;i++)
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
  if(fr->a[2]>=0)
    printf("Minimum in z direction, ");
  else
    printf("Maximum in z direction, ");
  if(fr->vertBoundsFound==1)
    {
      if ((float)(fr->fitVert[2]-fr->vertLBound[2])==(float)(fr->vertUBound[2]-fr->fitVert[2]))
        printf("z0 = %LE +/- %LE\n",fr->fitVert[2],fr->vertUBound[2]-fr->fitVert[2]);
      else
        printf("z0 = %LE + %LE - %LE\n",fr->fitVert[2],fr->vertUBound[2]-fr->fitVert[2],fr->fitVert[2]-fr->vertLBound[2]);
    }
  else
    printf("z0 = %LE\n",fr->fitVert[2]);
  
  printf("\nf(x0,y0,z0) = %LE\n",fr->vertVal); 
    
}

//generates the functional form of the fit function for plotting,
//which varies depending on the plotting mode (parameters may be fixed)
void plotForm3Par(const parameters * p, fit_results * fr, const plot_data * pd)
{
	//set up equation forms for plotting
	if(strcmp(p->plotMode,"1d")==0)
		{
			sprintf(fr->fitForm[0], "%LE*(x**2) + %LE*(%LE**2) + %LE*(%LE**2) + %LE*x*%LE + %LE*x*%LE + %LE*%LE*%LE + %LE*x + %LE*%LE + %LE*%LE + %LE",fr->a[0],fr->a[1],pd->fixedParVal[1],fr->a[2],pd->fixedParVal[2],fr->a[3],pd->fixedParVal[1],fr->a[4],pd->fixedParVal[2],fr->a[5],pd->fixedParVal[1],pd->fixedParVal[2],fr->a[6],fr->a[7],pd->fixedParVal[1],fr->a[8],pd->fixedParVal[2],fr->a[9]);
			sprintf(fr->fitForm[1], "%LE*(%LE**2) + %LE*(x**2) + %LE*(%LE**2) + %LE*x*%LE + %LE*%LE*%LE + %LE*x*%LE + %LE*x + %LE*%LE + %LE*%LE + %LE",fr->a[0],pd->fixedParVal[0],fr->a[1],fr->a[2],pd->fixedParVal[2],fr->a[3],pd->fixedParVal[0],fr->a[4],pd->fixedParVal[0],pd->fixedParVal[2],fr->a[5],pd->fixedParVal[2],fr->a[7],fr->a[6],pd->fixedParVal[0],fr->a[8],pd->fixedParVal[2],fr->a[9]);
			sprintf(fr->fitForm[2], "%LE*(%LE**2) + %LE*(%LE**2) + %LE*(x**2) + %LE*%LE*%LE + %LE*x*%LE + %LE*x*%LE + %LE*x + %LE*%LE + %LE*%LE + %LE",fr->a[0],pd->fixedParVal[0],fr->a[1],pd->fixedParVal[1],fr->a[2],fr->a[3],pd->fixedParVal[0],pd->fixedParVal[1],fr->a[4],pd->fixedParVal[0],fr->a[5],pd->fixedParVal[1],fr->a[8],fr->a[6],pd->fixedParVal[0],fr->a[7],pd->fixedParVal[1],fr->a[9]);
    }
  else if(strcmp(p->plotMode,"2d")==0)
		{
			sprintf(fr->fitForm[0], "%LE*(%LE**2) + %LE*(x**2) + %LE*(y**2) + %LE*%LE*x + %LE*%LE*y + %LE*x*y + %LE*%LE + %LE*x + %LE*y + %LE",fr->a[0],pd->fixedParVal[0],fr->a[1],fr->a[2],fr->a[3],pd->fixedParVal[0],fr->a[4],pd->fixedParVal[0],fr->a[5],fr->a[6],pd->fixedParVal[0],fr->a[7],fr->a[8],fr->a[9]);//x=par[1],y=par[2]
			sprintf(fr->fitForm[1], "%LE*(x**2) + %LE*(%LE**2) + %LE*(y**2) + %LE*x*%LE + %LE*x*y + %LE*%LE*y + %LE*x + %LE*%LE + %LE*y + %LE",fr->a[0],fr->a[1],pd->fixedParVal[1],fr->a[2],fr->a[3],pd->fixedParVal[1],fr->a[4],fr->a[5],pd->fixedParVal[1],fr->a[6],fr->a[7],pd->fixedParVal[1],fr->a[8],fr->a[9]);//x=par[0],y=par[2]
			sprintf(fr->fitForm[2], "%LE*(x**2) + %LE*(y**2) + %LE*(%LE**2) + %LE*x*y + %LE*x*%LE + %LE*y*%LE + %LE*x + %LE*y + %LE*%LE + %LE",fr->a[0],fr->a[1],fr->a[2],pd->fixedParVal[2],fr->a[3],fr->a[4],pd->fixedParVal[2],fr->a[5],pd->fixedParVal[2],fr->a[6],fr->a[7],fr->a[8],pd->fixedParVal[2],fr->a[9]);//x=par[0],y=par[1]
		}
}


//fit data to a paraboloid of the form
//f(x,y,z) = a1*x^2 + a2*y^2 + a3*z^2 + a4*x*y + a5*x*z + a6*y*z +a7*x + a8*y + a9*z + a10
void fit3Par(const parameters * p, const data * d, fit_results * fr, plot_data * pd, int print)
{
  //construct equations (n=3 specific case)
  int i,j;
  lin_eq_type linEq;
  linEq.dim=10;
  
  for(i=0;i<3;i++)//loop over free parameters
    for(j=i;j<3;j++)//loop over free parameters
      linEq.matrix[i][j]=d->xxpowsum[i][2][j][2];//top-left 3x3 entries
      
  linEq.matrix[0][3]=d->xxpowsum[0][3][1][1];
  linEq.matrix[0][4]=d->xxpowsum[0][3][2][1];
  linEq.matrix[0][5]=d->xxxpowsum[0][2][1][1][2][1];
  linEq.matrix[0][6]=d->xpowsum[0][3];
  linEq.matrix[0][7]=d->xxpowsum[0][2][1][1];
  linEq.matrix[0][8]=d->xxpowsum[0][2][2][1];
  linEq.matrix[0][9]=d->xpowsum[0][2];
  
  linEq.matrix[1][3]=d->xxpowsum[0][1][1][3];
  linEq.matrix[1][4]=d->xxxpowsum[0][1][1][2][2][1];
  linEq.matrix[1][5]=d->xxpowsum[1][3][2][1];
  linEq.matrix[1][6]=d->xxpowsum[0][1][1][2];
  linEq.matrix[1][7]=d->xpowsum[1][3];
  linEq.matrix[1][8]=d->xxpowsum[1][2][2][1];
  linEq.matrix[1][9]=d->xpowsum[1][2];
  
  linEq.matrix[2][3]=d->xxxpowsum[0][1][1][1][2][2];
  linEq.matrix[2][4]=d->xxpowsum[0][1][2][3];
  linEq.matrix[2][5]=d->xxpowsum[1][1][2][3];
  linEq.matrix[2][6]=d->xxpowsum[0][1][2][2];
  linEq.matrix[2][7]=d->xxpowsum[1][1][2][2];
  linEq.matrix[2][8]=d->xpowsum[2][3];
  linEq.matrix[2][9]=d->xpowsum[2][2];

  linEq.matrix[3][3]=d->xxpowsum[0][2][1][2];  
  linEq.matrix[3][4]=d->xxxpowsum[0][2][1][1][2][1];
  linEq.matrix[3][5]=d->xxxpowsum[0][1][1][2][2][1];
  linEq.matrix[3][6]=d->xxpowsum[0][2][1][1];
  linEq.matrix[3][7]=d->xxpowsum[0][1][1][2];
  linEq.matrix[3][8]=d->xxxpowsum[0][1][1][1][2][1];
  linEq.matrix[3][9]=d->xxpowsum[0][1][1][1];
  
  linEq.matrix[4][4]=d->xxpowsum[0][2][2][2];
  linEq.matrix[4][5]=d->xxxpowsum[0][1][1][1][2][2];
  linEq.matrix[4][6]=d->xxpowsum[0][2][2][1];
  linEq.matrix[4][7]=d->xxxpowsum[0][1][1][1][2][1];
  linEq.matrix[4][8]=d->xxpowsum[0][1][2][2];
  linEq.matrix[4][9]=d->xxpowsum[0][1][2][1];
  
  linEq.matrix[5][5]=d->xxpowsum[1][2][2][2];
  linEq.matrix[5][6]=d->xxxpowsum[0][1][1][1][2][1];
  linEq.matrix[5][7]=d->xxpowsum[1][2][2][1];
  linEq.matrix[5][8]=d->xxpowsum[1][1][2][2];
  linEq.matrix[5][9]=d->xxpowsum[1][1][2][1];
  
  linEq.matrix[6][6]=d->xpowsum[0][2];
  linEq.matrix[6][7]=d->xxpowsum[0][1][1][1];
  linEq.matrix[6][8]=d->xxpowsum[0][1][2][1];
  linEq.matrix[6][9]=d->xpowsum[0][1];
  
  linEq.matrix[7][7]=d->xpowsum[1][2];
  linEq.matrix[7][8]=d->xxpowsum[1][1][2][1];
  linEq.matrix[7][9]=d->xpowsum[1][1];
  
  linEq.matrix[8][8]=d->xpowsum[2][2];
  linEq.matrix[8][9]=d->xpowsum[2][1];
      
  linEq.matrix[9][9]=d->xpowsum[0][0];//bottom right entry
  
  //mirror the matrix (top right half mirrored to bottom left half)
  for(i=1;i<linEq.dim;i++)
    for(j=0;j<i;j++)
      linEq.matrix[i][j]=linEq.matrix[j][i];
  
  for(i=0;i<3;i++)
    linEq.vector[i]=d->mxpowsum[i][2];
  linEq.vector[3]=d->mxxpowsum[0][1][1][1];
  linEq.vector[4]=d->mxxpowsum[0][1][2][1];
  linEq.vector[5]=d->mxxpowsum[1][1][2][1];
  for(i=6;i<9;i++)
    linEq.vector[i]=d->mxpowsum[i-6][1];
  linEq.vector[9]=d->msum;
    
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
  fr->ndf=d->lines-11;
  for(i=0;i<d->lines;i++)//loop over data points for chisq
    {
      f=fr->a[0]*d->x[0][i]*d->x[0][i] + fr->a[1]*d->x[1][i]*d->x[1][i] + fr->a[2]*d->x[2][i]*d->x[2][i] + fr->a[3]*d->x[0][i]*d->x[1][i] + fr->a[4]*d->x[0][i]*d->x[2][i] + fr->a[5]*d->x[1][i]*d->x[2][i] + fr->a[6]*d->x[0][i] + fr->a[7]*d->x[1][i] + fr->a[8]*d->x[2][i] + fr->a[9];
      fr->chisq+=(d->x[3][i] - f)*(d->x[3][i] - f)/(d->x[3+1][i]*d->x[3+1][i]);
    }
  //Calculate covariances and uncertainties, see J. Wolberg 
  //'Data Analysis Using the Method of Least Squares' sec 2.5
  for(i=0;i<linEq.dim;i++)
    for(j=0;j<linEq.dim;j++)
      fr->covar[i][j]=linEq.inv_matrix[i][j]*(fr->chisq/fr->ndf);
  for(i=0;i<linEq.dim;i++)
    fr->aerr[i]=(long double)sqrt((double)(fr->covar[i][i]));
    
  //now that the fit is performed, use the fit parameters (and the derivative of the fitting function) to find the minimum
  linEq.dim=3;
  linEq.matrix[0][0]=2*fr->a[0];
  linEq.matrix[0][1]=fr->a[3];
  linEq.matrix[0][2]=fr->a[4];
  linEq.matrix[1][1]=2*fr->a[1];
  linEq.matrix[1][2]=fr->a[5];
  linEq.matrix[2][2]=2*fr->a[2];
  //mirror the matrix (top right half mirrored to bottom left half)
  for(i=1;i<linEq.dim;i++)
    for(j=0;j<i;j++)
      linEq.matrix[i][j]=linEq.matrix[j][i];
      
  linEq.vector[0]=-1*fr->a[6];
  linEq.vector[1]=-1*fr->a[7];
  linEq.vector[2]=-1*fr->a[8];
  
  //solve system of equations and assign values
  if(!(solve_lin_eq(&linEq)==1))
    {
      printf("ERROR: Could not determine paraboloid center point.\n");
      exit(-1);
    }
  
  for(i=0;i<linEq.dim;i++)
    fr->fitVert[i]=linEq.solution[i];   
  
	//find the value of the fit function at the vertex
	fr->vertVal=eval3Par(fr->fitVert[0],fr->fitVert[1],fr->fitVert[2],fr);
	
	if(strcmp(p->dataType,"chisq")==0)
		fit3ParChisqConf(p,fr);//generate confidence interval bounds for chisq data
	
	//print results
  if(print==1)
		print3Par(d,p,fr);
		
	if((p->plotData==1)&&(p->verbose<1))
		{
			preparePlotData(d,p,fr,pd);
			plotForm3Par(p,fr,pd);
			plotData(p,fr,pd);
		}
  
}
