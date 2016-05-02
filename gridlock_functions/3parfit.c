//fit data to a paraboloid of the form
//f(x,y,z) = a1*x^2 + a2*y^2 + a3*z^2 + a4*x*y + a5*x*z + a6*y*z +a7*x + a8*y + a9*z + a10
void fit3Par(const data * d, fit_results * fr)
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
      exit(-1);
    }
  
  //save fit parameters  
  for(i=0;i<linEq.dim;i++)
    fr->a[i]=linEq.solution[i];
    
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
  
}
