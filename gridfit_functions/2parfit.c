//fit data to a paraboloid of the form
//f(x,y) = a1*x^2 + a2*y^2 + a3*x*y + a4*x + a5*y + a6
void fit2Par(const data * d, fit_results * fr)
{
  //construct equations (n=3 specific case)
  int i,j;
  lin_eq_type linEq;
  linEq.dim=6;
  
  for(i=0;i<2;i++)//loop over free parameters
    for(j=0;j<2;j++)//loop over free parameters
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
    
  //now that the fit is performed, use the fit parameters (and the derivative of the fitting function) to find the minimum
  linEq.dim=2;
  linEq.matrix[0][0]=2*fr->a[0];
  linEq.matrix[0][1]=fr->a[2];
  linEq.matrix[1][0]=fr->a[2];
  linEq.matrix[1][1]=2*fr->a[1];
      
  linEq.vector[0]=-1*fr->a[3];
  linEq.vector[1]=-1*fr->a[4];
  
  //solve system of equations and assign values
  if(!(solve_lin_eq(&linEq)==1))
    {
      printf("ERROR: Could not determine paraboloid center point.\n");
      exit(-1);
    }
  
  for(i=0;i<linEq.dim;i++)
    fr->fitVert[i]=linEq.solution[i];

}
