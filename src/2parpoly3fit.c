//forward declarations
void generateSums(data *,const parameters *);

//evaluates the fit function at the specified point
long double eval2ParPoly3(long double x,long double y, const fit_results * fr)
{
	return fr->a[0]*x*x*x + fr->a[1]*y*y*y
					+ fr->a[2]*x*x*y + fr->a[3]*x*y*y
					+ fr->a[4]*x*x + fr->a[5]*y*y
					+ fr->a[6]*x*y + fr->a[7]*x + fr->a[8]*y + fr->a[9];
}

//determine the minimum and maximum bounds in each variable by generating 2 
//polynomials, one in each variable
//done by taking the minimum value of the fit function available 
//for various values of the variable of interest (projecting the 
//minimum values on the variable axis)
//fixZero - 0: don't fix minimum to 0, 1: fix minimum in x to 0, 
//             2: fix minimum in y to 0, 3: fix minimum in x and y to 0 
void fit2ParPoly3ChisqConf(const data * d, const parameters * p, fit_results * fr, int fixZero)
{
  int i,j;
	
  //allocate fit structures
  parameters *svarp=(parameters*)calloc(1,sizeof(parameters));
  data *svard=(data*)calloc(1,sizeof(data));
  fit_results *svarfr=(fit_results*)calloc(1,sizeof(fit_results));
  plot_data *svarpd=(plot_data*)calloc(1,sizeof(plot_data));
  //setup fit
  svarp->numVar=1;
  svarp->ciDelta=p->ciDelta;
  strcpy(svarp->plotMode,"1d");
  if(strcmp(p->dataType,"chisq")==0)
    strcpy(svarp->dataType,"chisq");
  	
	long double val,a,b,c,sqrtval;
	for(i=0;i<2;i++)//variable #
		{
			svard->lines=0;
			for(j=0;j<=100;j++)//number of data points to compute
				{
					val=d->min_x[i] + (j/100.)*(d->max_x[i] - d->min_x[i]);
					if(i==0)
						{
							a=3.*fr->a[1];
							b=2.*fr->a[3]*val + 2.*fr->a[5];
							c=fr->a[2]*val*val + fr->a[6]*val + fr->a[8];
						}
					else
						{
							a=3.*fr->a[0];
							b=2.*fr->a[2]*val + 2.*fr->a[4];
							c=fr->a[3]*val*val + fr->a[6]*val + fr->a[7];
						}
					sqrtval=b*b - 4.*a*c;			
					if(sqrtval==0)
						{
							if(i==0)
								svard->x[1][svard->lines]=eval2ParPoly3(val,-1.*b/(2.*a),fr);
							else
								svard->x[1][svard->lines]=eval2ParPoly3(-1.*b/(2.*a),val,fr);
						}
					else if(sqrtval>0)
						{
							if(i==0)
								{
									if(eval2ParPoly3(val,(-1.*b - sqrt(sqrtval))/(2.*a),fr)>eval2ParPoly3(val,(-1.*b + sqrt(sqrtval))/(2.*a),fr))
										svard->x[1][svard->lines]=eval2ParPoly3(val,(-1.*b + sqrt(sqrtval))/(2.*a),fr);
									else
										svard->x[1][svard->lines]=eval2ParPoly3(val,(-1.*b - sqrt(sqrtval))/(2.*a),fr);
								}
							else
								{
									if(eval2ParPoly3((-1.*b - sqrt(sqrtval))/(2.*a),val,fr)>eval2ParPoly3((-1.*b + sqrt(sqrtval))/(2.*a),val,fr))
										svard->x[1][svard->lines]=eval2ParPoly3((-1.*b + sqrt(sqrtval))/(2.*a),val,fr);
									else
										svard->x[1][svard->lines]=eval2ParPoly3((-1.*b - sqrt(sqrtval))/(2.*a),val,fr);
								}
						}//don't take any action (skip data point) if the roots are not real
					
					if(sqrtval>=0)
						{
							svard->x[0][svard->lines]=val;
							svard->x[2][svard->lines]=1.;//set weight
							svard->lines++;
						}
					
					
				}
				//fit and find critical points of this data
				generateSums(svard,svarp);
				fitPoly3(svarp,svard,svarfr,svarpd,0);//fit but don't print data

        //compute confidence interval with variable fixed to 0 if requested
        if(strcmp(p->dataType,"chisq")==0)
          if((fixZero==i+1)||(fixZero==3))
            fitPoly3ChisqConf(svarp,svarfr,0.);

				//save critical point corresponding to minimum
        if((fixZero==i+1)||(fixZero==3))
          fr->fitVert[i]=0.;
				else if(evalPoly3(svarfr->fitVert[0],svarfr)<evalPoly3(svarfr->fitVert[1],svarfr))
					fr->fitVert[i]=svarfr->fitVert[0];
				else
					fr->fitVert[i]=svarfr->fitVert[1];
				//save confidence bounds, if applicable
				if((strcmp(p->dataType,"chisq")==0)&&(svarfr->vertBoundsFound==1))
					{
						fr->vertBoundsFound=1;
						fr->vertLBound[i]=svarfr->vertLBound[0];
						fr->vertUBound[i]=svarfr->vertUBound[0];
					}
				else
					fr->vertBoundsFound=0;
				
					
			}
	
	//free fit structures
	free(svarp);
	free(svard);
	free(svarfr);
	free(svarpd);
}

void printFitVertex2ParPoly3(const data * d, const parameters * p, const fit_results * fr)
{
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

//prints fit data
void print2ParPoly3(const data * d, const parameters * p, const fit_results * fr)
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
  printf("Fit parameter uncertainties reported at 1-sigma.\n");
  printf("Fit function: f(x,y,z) = a1*x^3 + a2*y^3 + a3*x^2*y\n                       + a4*x*y^2 + a5*x^2 + a6*y^2\n                       + a7*x*y + a8*x + a9*y + a10\n\n");
  //printf("Best chisq (fit): %0.3Lf\nBest chisq/NDF (fit): %0.3Lf\n\n",fr->chisq,fr->chisq/fr->ndf);
  printf("Coefficients from fit: a1 = %LE +/- %LE\n",fr->a[0],fr->aerr[0]);
  for(i=1;i<10;i++)
    printf("                       a%i = %LE +/- %LE\n",i+1,fr->a[i],fr->aerr[i]);
  printf("\n");
  
  //print local minimum and confidence bounds (if necessary)
	if((strcmp(p->dataType,"chisq")==0)&&(fr->vertBoundsFound==1))
		{
      printf("Local minimum with %s confidence interval:\n",p->ciSigmaDesc);
			if((float)(fr->vertUBound[0]-fr->fitVert[0])==(float)(fr->fitVert[0]-fr->vertLBound[0]))
				printf("x = %LE +/- %LE\n",fr->fitVert[0],fr->vertUBound[0]-fr->fitVert[0]);
			else
				printf("x = %LE + %LE - %LE\n",fr->fitVert[0],fr->vertUBound[0]-fr->fitVert[0],fr->fitVert[0]-fr->vertLBound[0]);
			if((float)(fr->vertUBound[1]-fr->fitVert[1])==(float)(fr->fitVert[1]-fr->vertLBound[1]))
				printf("y = %LE +/- %LE\n",fr->fitVert[1],fr->vertUBound[1]-fr->fitVert[1]);
			else
				printf("y = %LE + %LE - %LE\n",fr->fitVert[1],fr->vertUBound[1]-fr->fitVert[1],fr->fitVert[1]-fr->vertLBound[1]);
		}
	else
    {
      printf("Local minimum:\n");
      printf("x = %LE, y = %LE\n",fr->fitVert[0],fr->fitVert[1]);
    }
		
    
}

void plotForm2ParPoly3(const parameters * p, fit_results * fr, const plot_data * pd)
{
	//set up equation forms for plotting
	if(strcmp(p->plotMode,"1d")==0)
		{
			sprintf(fr->fitForm[0], "%.12LE*(x**3) + %.12LE*(%.12LE**3) + %.12LE*(x**2)*%.12LE + %.12LE*x*(%.12LE**2) + %.12LE*(x**2) + %.12LE*(%.12LE**2) + %.12LE*x*%.12LE + %.12LE*x + %.12LE*%.12LE + %.12LE",fr->a[0],fr->a[1],pd->fixedParVal[1],fr->a[2],pd->fixedParVal[1],fr->a[3],pd->fixedParVal[1],fr->a[4],fr->a[5],pd->fixedParVal[1],fr->a[6],pd->fixedParVal[1],fr->a[7],fr->a[8],pd->fixedParVal[1],fr->a[9]);//y fixed
			sprintf(fr->fitForm[1], "%.12LE*(x**3) + %.12LE*(%.12LE**3) + %.12LE*(x**2)*%.12LE + %.12LE*x*(%.12LE**2) + %.12LE*(x**2) + %.12LE*(%.12LE**2) + %.12LE*x*%.12LE + %.12LE*x + %.12LE*%.12LE + %.12LE",fr->a[1],fr->a[0],pd->fixedParVal[0],fr->a[3],pd->fixedParVal[0],fr->a[2],pd->fixedParVal[0],fr->a[5],fr->a[4],pd->fixedParVal[0],fr->a[6],pd->fixedParVal[0],fr->a[8],fr->a[7],pd->fixedParVal[0],fr->a[9]);//x fixed
    }
  else if(strcmp(p->plotMode,"2d")==0)
    {
      sprintf(fr->fitForm[0], "%.12LE*(x**3) + %.12LE*(y**3) + %.12LE*(x**2)*y + %.12LE*(y**2)*x + %.12LE*(x**2) + %.12LE*(y**2) + %.12LE*x*y + %.12LE*x + %.12LE*y + %.12LE",fr->a[0],fr->a[1],fr->a[2],fr->a[3],fr->a[4],fr->a[5],fr->a[6],fr->a[7],fr->a[8],fr->a[9]);
    }
}

//fit data to a paraboloid of the form
//f(x,y) = a1*x^3 + a2*y^3 + a3*x^2*y + a4*x*y^2 + a5*x^2 + a6*y^2 +a7*x*y + a8*x + a9*y + a10
void fit2ParPoly3(const parameters * p, const data * d, fit_results * fr, plot_data * pd, int print)
{
  //construct equations
  int i,j;
  lin_eq_type linEq;
  linEq.dim=10;
  
  linEq.matrix[0][0]=d->xxpowsum[0][3][0][3];//x^6
  linEq.matrix[0][1]=d->xxpowsum[0][3][1][3];//x^3y^3
  linEq.matrix[0][2]=d->xxpowsum[0][5][1][1];
  linEq.matrix[0][3]=d->xxpowsum[0][4][1][2];
  linEq.matrix[0][4]=d->xpowsum[0][5];//x^5
  linEq.matrix[0][5]=d->xxpowsum[0][3][1][2];
  linEq.matrix[0][6]=d->xxpowsum[0][4][1][1];
  linEq.matrix[0][7]=d->xpowsum[0][4];
  linEq.matrix[0][8]=d->xxpowsum[0][3][1][1];
  linEq.matrix[0][9]=d->xpowsum[0][3];
  
  linEq.matrix[1][1]=d->xxpowsum[1][3][1][3];//y^6
  linEq.matrix[1][2]=d->xxpowsum[0][2][1][4];
  linEq.matrix[1][3]=d->xxpowsum[0][1][1][5];
  linEq.matrix[1][4]=d->xxpowsum[0][2][1][3];
  linEq.matrix[1][5]=d->xpowsum[1][5];
  linEq.matrix[1][6]=d->xxpowsum[0][1][1][4];
  linEq.matrix[1][7]=d->xxpowsum[0][1][1][3];
  linEq.matrix[1][8]=d->xpowsum[1][4];
  linEq.matrix[1][9]=d->xpowsum[1][3];
  
  linEq.matrix[2][2]=d->xxpowsum[0][4][1][2];
  linEq.matrix[2][3]=d->xxpowsum[0][3][1][3];
  linEq.matrix[2][4]=d->xxpowsum[0][4][1][1];
  linEq.matrix[2][5]=d->xxpowsum[0][2][1][3];
  linEq.matrix[2][6]=d->xxpowsum[0][3][1][2];
  linEq.matrix[2][7]=d->xxpowsum[0][3][1][1];
  linEq.matrix[2][8]=d->xxpowsum[0][2][1][2];
  linEq.matrix[2][9]=d->xxpowsum[0][2][1][1];

  linEq.matrix[3][3]=d->xxpowsum[0][2][1][4];  
  linEq.matrix[3][4]=d->xxpowsum[0][3][1][2];
  linEq.matrix[3][5]=d->xxpowsum[0][1][1][4];
  linEq.matrix[3][6]=d->xxpowsum[0][2][1][3];
  linEq.matrix[3][7]=d->xxpowsum[0][2][1][2];
  linEq.matrix[3][8]=d->xxpowsum[0][1][1][3];
  linEq.matrix[3][9]=d->xxpowsum[0][1][1][2];
  
  linEq.matrix[4][4]=d->xpowsum[0][4];
  linEq.matrix[4][5]=d->xxpowsum[0][2][1][2];
  linEq.matrix[4][6]=d->xxpowsum[0][3][1][1];
  linEq.matrix[4][7]=d->xpowsum[0][3];
  linEq.matrix[4][8]=d->xxpowsum[0][2][1][1];
  linEq.matrix[4][9]=d->xpowsum[0][2];
  
  linEq.matrix[5][5]=d->xpowsum[1][4];
  linEq.matrix[5][6]=d->xxpowsum[0][1][1][3];
  linEq.matrix[5][7]=d->xxpowsum[0][1][1][2];
  linEq.matrix[5][8]=d->xpowsum[1][3];
  linEq.matrix[5][9]=d->xpowsum[1][2];
  
  linEq.matrix[6][6]=d->xxpowsum[0][2][1][2];
  linEq.matrix[6][7]=d->xxpowsum[0][2][1][1];
  linEq.matrix[6][8]=d->xxpowsum[0][1][1][2];
  linEq.matrix[6][9]=d->xxpowsum[0][1][1][1];
  
  linEq.matrix[7][7]=d->xpowsum[0][2];
  linEq.matrix[7][8]=d->xxpowsum[0][1][1][1];
  linEq.matrix[7][9]=d->xpowsum[0][1];
  
  linEq.matrix[8][8]=d->xpowsum[1][2];
  linEq.matrix[8][9]=d->xpowsum[1][1];
      
  linEq.matrix[9][9]=d->xpowsum[0][0];//bottom right entry
  
  //mirror the matrix (top right half mirrored to bottom left half)
  for(i=1;i<linEq.dim;i++)
    for(j=0;j<i;j++)
      linEq.matrix[i][j]=linEq.matrix[j][i];
  
  linEq.vector[0]=d->mxpowsum[0][3];
  linEq.vector[1]=d->mxpowsum[1][3];
  linEq.vector[2]=d->mxxpowsum[0][2][1][1];
  linEq.vector[3]=d->mxxpowsum[0][1][1][2];
  linEq.vector[4]=d->mxpowsum[0][2];
  linEq.vector[5]=d->mxpowsum[1][2];
  linEq.vector[6]=d->mxxpowsum[0][1][1][1];
  linEq.vector[7]=d->mxpowsum[0][1];
  linEq.vector[8]=d->mxpowsum[1][1];
  linEq.vector[9]=d->msum;

  /*printf("Matrix:\n");
  for(i=0;i<linEq.dim;i++)
    {
      for(j=0;j<linEq.dim;j++)
        printf(" %LE ",linEq.matrix[i][j]);
      printf("\n");
    }
  printf("Vector:\n");
  for(i=0;i<linEq.dim;i++)
    {
      printf("%LE\n",linEq.vector[i]);
    }*/

	//solve system of equations and assign values
	if(!(solve_lin_eq(&linEq)==1))
		{
			printf("ERROR: Could not determine fit parameters (2parpoly3).\n");
			printf("Perhaps there are not enough data points to perform a fit?\n");
      printf("Otherwise you can also try adjusting the fit range using the UPPER_LIMITS and LOWER_LIMITS options.\n");
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
      f=fr->a[0]*d->x[0][i]*d->x[0][i]*d->x[0][i] + fr->a[1]*d->x[1][i]*d->x[1][i]*d->x[1][i]
       + fr->a[2]*d->x[0][i]*d->x[0][i]*d->x[1][i] + fr->a[3]*d->x[0][i]*d->x[1][i]*d->x[1][i]
        + fr->a[4]*d->x[0][i]*d->x[0][i] + fr->a[5]*d->x[1][i]*d->x[1][i]
         + fr->a[6]*d->x[0][i]*d->x[1][i] + fr->a[7]*d->x[0][i] + fr->a[8]*d->x[1][i] + fr->a[9];
      fr->chisq+=(d->x[2][i] - f)*(d->x[2][i] - f)/(d->x[2+1][i]*d->x[2+1][i]);
    }
  
  /*printf("Inverse matrix:\n");
  for(i=0;i<linEq.dim;i++)
    {
      for(j=0;j<linEq.dim;j++)
        printf(" %LE ",linEq.inv_matrix[i][j]);
      printf("\n");
    }
  printf("chisq: %Lf\n",fr->chisq);*/

  //Calculate covariances and uncertainties, see J. Wolberg 
  //'Data Analysis Using the Method of Least Squares' sec 2.5
  for(i=0;i<linEq.dim;i++)
    for(j=0;j<linEq.dim;j++)
      fr->covar[i][j]=linEq.inv_matrix[i][j]*(fr->chisq/fr->ndf);
  for(i=0;i<linEq.dim;i++)
    fr->aerr[i]=(long double)sqrt((double)(fr->covar[i][i]));
	
	
	fit2ParPoly3ChisqConf(d,p,fr,0);
	

	//print results
  if(print==1)
    print2ParPoly3(d,p,fr);
  
  //check zero bounds
  if(strcmp(p->dataType,"chisq")==0)
    {
      if((fr->fitVert[0]<0.) || (p->forceZeroX==1))
        {
          if((fr->fitVert[1]>=0.) || (p->forceZeroX==1))
            {
              //make a copy of the fit results to work on
              fit_results *temp1=(fit_results*)calloc(1,sizeof(fit_results));
              memcpy(temp1,fr,sizeof(fit_results));
              //compute fit vertex assuming x is fixed to 0 (from 1st derivative condition)
              temp1->fitVert[0]=0.;
              long double yVertCandidate = (-2.*temp1->a[5] + sqrt((2.*temp1->a[5])*(2.*temp1->a[5]) - 12.*temp1->a[1]*temp1->a[8]))/(6.*temp1->a[1]);
              //check concavity (from 2nd derivative condition)
              if(yVertCandidate > (-2.*temp1->a[5]/(6.*temp1->a[1])))
                temp1->fitVert[1] = yVertCandidate;
              else
                temp1->fitVert[1] = (-2.*temp1->a[5] - sqrt((2.*temp1->a[5])*(2.*temp1->a[5]) - 12.*temp1->a[1]*temp1->a[8]))/(6.*temp1->a[1]);
              temp1->vertVal=eval2ParPoly3(temp1->fitVert[0],temp1->fitVert[1],temp1);
              //fit confidence interval to get bound for x
              fit2ParPoly3ChisqConf(d,p,temp1,1);
              temp1->vertLBound[0]=-1.*temp1->vertUBound[0]; //mirror bounds in x
              //determine confidence in y by fitting the 3rd order polynomial at x=0
              fit_results *temp2=(fit_results*)calloc(1,sizeof(fit_results));
              //map fit parameters into a 3rd order polynomial
              temp2->a[0]=temp1->a[1];
              temp2->a[1]=temp1->a[5];
              temp2->a[2]=temp1->a[8];
              temp2->a[3]=temp1->a[9];
              //fit confidence interval of the 3rd order polynomial to get bound for y (using ciDelta assuming 2 variables)
              fitPoly3ChisqConf(p,temp2,temp1->fitVert[1]);
              temp1->vertLBound[1]=temp2->vertLBound[0];
              temp1->vertUBound[1]=temp2->vertUBound[0];
              free(temp2);
              if(print==1)
                {
                  printf("\nAssuming minimum at zero for x,\n");
                  printFitVertex2ParPoly3(d,p,temp1);
                }
              free(temp1);
            }
          
        }
      else if((fr->fitVert[1]<0.) || (p->forceZeroY==1)) //implies fitVert[0]>=0.
        {
          //make a copy of the fit results to work on
          fit_results *temp1=(fit_results*)calloc(1,sizeof(fit_results));
          memcpy(temp1,fr,sizeof(fit_results));
          //compute fit vertex assuming y is fixed to 0 (from 1st derivative condition)
          temp1->fitVert[1]=0.;
          long double xVertCandidate = (-2.*temp1->a[4] + sqrt((2.*temp1->a[4])*(2.*temp1->a[4]) - 12.*temp1->a[0]*temp1->a[7]))/(6.*temp1->a[0]);
          if(xVertCandidate > (-2.*temp1->a[4]/(6.*temp1->a[0])))
            temp1->fitVert[0] = xVertCandidate;
          else
            temp1->fitVert[0] = (-2.*temp1->a[4] - sqrt((2.*temp1->a[4])*(2.*temp1->a[4]) - 12.*temp1->a[0]*temp1->a[7]))/(6.*temp1->a[0]);
          temp1->vertVal=eval2ParPoly3(temp1->fitVert[0],temp1->fitVert[1],temp1);
          //fit confidence interval to get bound for y
          fit2ParPoly3ChisqConf(d,p,temp1,2);
          temp1->vertLBound[1]=-1.*temp1->vertUBound[1]; //mirror bounds in y
          //determine confidence in x by fitting the 3rd order polynomial at y=0
          fit_results *temp2=(fit_results*)calloc(1,sizeof(fit_results));
          //map fit parameters into a 3rd order polynomial
          temp2->a[0]=temp1->a[0];
          temp2->a[1]=temp1->a[4];
          temp2->a[2]=temp1->a[7];
          temp2->a[3]=temp1->a[9];
          //fit confidence interval of the 3rd order polynomial to get bound for x (using ciDelta assuming 2 variables)
          fitPoly3ChisqConf(p,temp2,temp1->fitVert[0]);
          temp1->vertLBound[0]=temp2->vertLBound[0];
          temp1->vertUBound[0]=temp2->vertUBound[0];
          free(temp2);
          if(print==1)
            {
              printf("\nAssuming minimum at zero for y,\n");
              printFitVertex2ParPoly3(d,p,temp1);
            }
          free(temp1);
        }
    }
		

	
	if((p->plotData==1)&&(p->verbose<1))
		{
			preparePlotData(d,p,fr,pd);
			plotForm2ParPoly3(p,fr,pd);
			plotData(p,fr,pd);
		}

}
