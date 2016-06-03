//fit data to a line of the form
//f(x) = a1*x + a2
//assuming errors in both x and y and perpendicular to the slope
//as this problem cannot be linearized it will be fit iteravely
void fitLinSlopeErrors(const data * d, fit_results * fr)
{
	//first, fit data to line assuming only y-errors
	//this is done to get an initial guess for the slope and intercept
  //construct equations
  int i,j,k;
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
  
  //save initial fit parameters  
  for(i=0;i<linEq.dim;i++)
    fr->a[i]=linEq.solution[i];
  
  
  
  //set up grid to iterate over
  //this is a 2D grid in the slope and intercept
  parameters *gridPar=(parameters*)calloc(1,sizeof(parameters));
  data *gridPts=(data*)calloc(1,sizeof(data));
  fit_results *gridRes=(fit_results*)calloc(1,sizeof(fit_results));
  gridPar->numVar=2;
  //grid search parameters
  int gridDim=10;
  double frac=0.1;//full fractional range to vary grid parameters by from initial guess
  //chisq parameters
  long double f;
	fr->ndf=d->lines-3;
  
  printf("initial slope	intercept: %Lf %Lf\n",fr->a[0],fr->a[1]);
  printf("slope	intercept	chisq\n");
  int pt=0;
  for(i=0;i<gridDim;i++)
  	for(j=0;j<gridDim;j++)
  	{
  		gridPts->x[0][pt]= fr->a[0] + (i/(double)gridDim)*frac*fr->a[0] - 0.5*frac*fr->a[0];
  		gridPts->x[1][pt]= fr->a[1] + (j/(double)gridDim)*frac*fr->a[1] - 0.5*frac*fr->a[1];
  		
  		//get chisq (NEED TO IMPLEMENT PROPER MEASURE HERE USING PERPENDICULAR LINE)
  		fr->chisq=0;
			for(k=0;k<d->lines;k++)//loop over data points for chisq
				{
				  f=gridPts->x[0][pt]*d->x[0][k] + gridPts->x[1][pt];
				  fr->chisq+=(d->x[1][k] - f)*(d->x[1][k] - f)
				              /(d->x[1+1][k]*d->x[1+1][k]);
				}
  		
  		gridPts->x[2][pt]=fr->chisq/fr->ndf;
  		printf("%Lf %Lf %Lf\n",gridPts->x[0][pt],gridPts->x[1][pt],gridPts->x[2][pt]);
  		pt++;
  	}
  generateSums(gridPts,gridPar);
  fit2Par(gridPts,gridRes);
  fit2ParChisqConf(gridRes);//generate confidence interval bounds for chisq data
  
	//save parameters from grid fit
	for(i=0;i<linEq.dim;i++)
  	{
			fr->a[i]=gridRes->fitVert[i];
			fr->aerr[i]=gridRes->vertUBound[i];
		}
	
	//now that the fit is performed, use the fit parameters to find the intercept(s)
  fr->fitVert[0]=-1.0*fr->a[1]/fr->a[0];//x-intercept
  fr->fitVert[1]=fr->a[1];//y-intercept
  
  //free structures
  free(gridPar);
  free(gridPts);
  free(gridRes);
}
