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
    
  //now that the fit is performed, use the fit parameters (and the derivative of the fitting function) to find the minimum
  fr->fitVert[0]=-1.0*fr->a[1]/(2.*fr->a[0]);

}
