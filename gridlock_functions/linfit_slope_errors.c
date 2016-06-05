//fit data to a line of the form
//f(x) = a1*x + a2
//assuming errors in both x and y and perpendicular to the slope
//as this problem cannot be linearized it will be fit iteravely
void fitLinSlopeErrors(const data * d, fit_results * fr)
{

  int i,j,k;

	//first, fit data to line assuming only y-errors
	//this is done to get an initial guess for the slope and intercept
  //construct equations
  fitLin(d,fr);
  
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
	long double r2,ls,li;//residualsq, line slope, line intercept
  long double ix,iy;//intersection x,intersection y
  
  for(i=0;i<=gridDim;i++)
  	for(j=0;j<=gridDim;j++)
    	{
    		gridPts->x[0][gridPts->lines]= fr->a[0] + (i/(double)gridDim)*frac*fr->a[0] - 0.5*frac*fr->a[0];
    		gridPts->x[1][gridPts->lines]= fr->a[1] + (j/(double)gridDim)*frac*fr->a[1] - 0.5*frac*fr->a[1];
    		
    		//get chisq
    		fr->chisq=0;
			  for(k=0;k<d->lines;k++)//loop over data points for chisq
				  {
				    //get residual which is perpendicular to fit line
				    ls=-1./gridPts->x[0][gridPts->lines];
				    li=d->x[1][k] - ls*d->x[0][k];
				    ix=(gridPts->x[1][gridPts->lines]-li)/(ls-gridPts->x[0][gridPts->lines]);
				    iy=ls*ix + li;
				    r2=(d->x[0][k]-ix)*(d->x[0][k]-ix) + (d->x[1][k]-iy)*(d->x[1][k]-iy);
				    fr->chisq+=r2/(d->x[1+1][k]*d->x[1+1][k]);
				  }
    		
    		gridPts->x[2][gridPts->lines]=fr->chisq/fr->ndf;
    		gridPts->x[3][gridPts->lines]=1.0;//need to assign weights
    		//printf("%Lf %Lf %Lf\n",gridPts->x[0][gridPts->lines],gridPts->x[1][gridPts->lines],gridPts->x[2][gridPts->lines]);
    		gridPts->lines++;
    	}

  generateSums(gridPts,gridPar);
  fit2Par(gridPts,gridRes);
  fit2ParChisqConf(gridRes);//generate confidence interval bounds for chisq data
  
	//save parameters from grid fit
	for(i=0;i<2;i++)
  	{
			fr->a[i]=gridRes->fitVert[i];
			fr->aerr[i]=gridRes->vertUBound[i];
		}
	
	//get chisq
	fr->chisq=0;
  for(i=0;i<d->lines;i++)//loop over data points for chisq
	  {
	    //get residual which is perpendicular to fit line
	    ls=-1./fr->a[0];
	    li=d->x[1][i] - ls*d->x[0][i];
	    ix=(fr->a[1]-li)/(ls-fr->a[0]);
	    iy=ls*ix + li;
	    r2=(d->x[0][i]-ix)*(d->x[0][i]-ix) + (d->x[1][i]-iy)*(d->x[1][i]-iy);
	    fr->chisq+=r2/(d->x[1+1][i]*d->x[1+1][i]);
	  }
	
	//now that the fit is performed, use the fit parameters to find the intercept(s)
  fr->fitVert[0]=-1.0*fr->a[1]/fr->a[0];//x-intercept
  fr->fitVert[1]=fr->a[1];//y-intercept
  
  //free structures
  free(gridPar);
  free(gridPts);
  free(gridRes);
}
