//prints fit data
void printResults(const data * d, const parameters * p, const fit_results * fr)
{

  int i;

  //simplified data printing depending on verbosity setting
  if(p->verbose==1)
    {
      for(i=0;i<p->numVar;i++)
        printf("%LE ",fr->fitVert[i]);
      printf("\n");
      return;
    }
  
  if(strcmp(p->fitType,"lin")==0)
    {
      printf("\nFIT RESULTS\n-----------\n");
      printf("Uncertainties reported at 1-sigma.\n");
      printf("Fit function: f(x,y) = a1*x + a2\n\n");
      printf("Best chisq (fit): %0.3Lf\nBest chisq/NDF (fit): %0.3Lf\n\n",fr->chisq,fr->chisq/fr->ndf);
      printf("Coefficients from fit: a1 = %LE +/- %LE\n",fr->a[0],fr->aerr[0]);
      printf("                       a2 = %LE +/- %LE\n",fr->a[1],fr->aerr[1]);
      printf("\n");
      
      printf("x-intercept = %LE\n",fr->fitVert[0]);
      printf("y-intercept = %LE\n",fr->fitVert[1]);

    }

  if(strcmp(p->fitType,"par1")==0)
    {
      printf("\nFIT RESULTS\n-----------\n");
      printf("Uncertainties reported at 1-sigma.\n");
      printf("Fit function: f(x,y) = a1*x^2 + a2*x + a3\n\n");
      printf("Best chisq (fit): %0.3Lf\nBest chisq/NDF (fit): %0.3Lf\n\n",fr->chisq,fr->chisq/fr->ndf);
      printf("Coefficients from fit: a1 = %LE +/- %LE\n",fr->a[0],fr->aerr[0]);
      for(i=1;i<3;i++)
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
      
      printf("\nf(x0) = %LE\n",fr->vertVal); 

    }

  if(strcmp(p->fitType,"par2")==0)
    {
    
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

  if(strcmp(p->fitType,"par3")==0)
    {
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
  
}
