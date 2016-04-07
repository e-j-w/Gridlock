//prints fit data
void printResults(const data * d, const parameters * p, const fit_results * fr)
{

  int i;

  if(p->numVar==2)
    {
      printf("\nFIT RESULTS\n-----------\n");
      printf("Fit function: f(x,y) = a1*x^2 + a2*y^2 + a3*x*y\n                     + a4*x + a5*y + a6\n\n");
      printf("Coefficients from fit: a1 = %LE\n",fr->a[0]);
      for(i=1;i<6;i++)
        printf("                       a%i = %LE\n",i+1,fr->a[i]);
      printf("\n");
      
      if(fr->a[0]>=0)
        printf("Minimum in x direction, ");
      else
        printf("Maximum in x direction, ");
      printf("x0 = %LE\n",fr->fitVert[0]);
      if(fr->a[1]>=0)
        printf("Minimum in y direction, ");
      else
        printf("Maximum in y direction, ");
      printf("y0 = %LE\n",fr->fitVert[1]);
      
      long double fitVal=fr->a[0]*fr->fitVert[0]*fr->fitVert[0] + fr->a[1]*fr->fitVert[1]*fr->fitVert[1] + fr->a[2]*fr->fitVert[0]*fr->fitVert[1] + fr->a[3]*fr->fitVert[0] + fr->a[4]*fr->fitVert[1] + fr->a[5];
      
      printf("\nf(x0,y0) = %LE\n",fitVal); 
      
      //get chisq value
      long double f;
      long double chisq=0;
      int ndf=d->lines-7;
      for(i=0;i<d->lines;i++)//loop over data points
        {
          f=fr->a[0]*d->x[0][i]*d->x[0][i] + fr->a[1]*d->x[1][i]*d->x[1][i] + fr->a[2]*d->x[0][i]*d->x[1][i] + fr->a[3]*d->x[0][i] + fr->a[4]*d->x[1][i] + fr->a[5];
          chisq+=(d->x[p->numVar][i] - f)*(d->x[p->numVar][i] - f);
        }
      printf("\nchisq: %0.3Lf\nchisq/NDF: %0.3Lf\n",chisq,chisq/ndf);
    }

  if(p->numVar==3)
    {
      printf("\nFIT RESULTS\n-----------\n");
      printf("Fit function: f(x,y,z) = a1*x^2 + a2*y^2 + a3*z^2\n                       + a4*x*y + a5*x*z + a6*y*z\n                       + a7*x + a8*y + a9*z + a10\n\n");
      printf("Coefficients from fit: a1 = %LE\n",fr->a[0]);
      for(i=1;i<10;i++)
        printf("                       a%i = %LE\n",i+1,fr->a[i]);
      printf("\n");
      
      if(fr->a[0]>=0)
        printf("Minimum in x direction, ");
      else
        printf("Maximum in x direction, ");
      printf("x0 = %LE\n",fr->fitVert[0]);
      if(fr->a[1]>=0)
        printf("Minimum in y direction, ");
      else
        printf("Maximum in y direction, ");
      printf("y0 = %LE\n",fr->fitVert[1]);
      if(fr->a[2]>=0)
        printf("Minimum in z direction, ");
      else
        printf("Maximum in z direction, ");  
      printf("z0 = %LE\n",fr->fitVert[2]);
      
      long double fitVal=fr->a[0]*fr->fitVert[0]*fr->fitVert[0] + fr->a[1]*fr->fitVert[1]*fr->fitVert[1] + fr->a[2]*fr->fitVert[2]*fr->fitVert[2] + fr->a[3]*fr->fitVert[0]*fr->fitVert[1] + fr->a[4]*fr->fitVert[0]*fr->fitVert[2] + fr->a[5]*fr->fitVert[1]*fr->fitVert[2] + fr->a[6]*fr->fitVert[0] + fr->a[7]*fr->fitVert[1] + fr->a[8]*fr->fitVert[2] + fr->a[9];
      
      printf("\nf(x0,y0,z0) = %LE\n",fitVal);
      
      //get chisq value
      long double f;
      long double chisq=0;
      int ndf=d->lines-11;
      for(i=0;i<d->lines;i++)//loop over data points
        {
          f=fr->a[0]*d->x[0][i]*d->x[0][i] + fr->a[1]*d->x[1][i]*d->x[1][i] + fr->a[2]*d->x[2][i]*d->x[2][i] + fr->a[3]*d->x[0][i]*d->x[1][i] + fr->a[4]*d->x[0][i]*d->x[2][i] + fr->a[5]*d->x[1][i]*d->x[2][i] + fr->a[6]*d->x[0][i] + fr->a[7]*d->x[1][i] + fr->a[8]*d->x[2][i] + fr->a[9];
          chisq+=(d->x[p->numVar][i] - f)*(d->x[p->numVar][i] - f);
        }
        
      printf("\nchisq: %0.3Lf\nchisq/NDF: %0.3Lf\n",chisq,chisq/ndf);
    } 
  
}
