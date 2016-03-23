//fit data to a paraboloid of the form
//f(x,y) = a1*x^2 + a2*y^2 + a3*x*y + a4*x + a5*y + a6
void fit2Par()
{
  //construct equations (n=3 specific case)
  linEq.dim=6;
  
  for(i=0;i<2;i++)//loop over free parameters
    for(j=0;j<2;j++)//loop over free parameters
      linEq.matrix[i][j]=xxpowsum[i][2][j][2];//top-left 2x2 entries
      
  linEq.matrix[0][2]=xxpowsum[0][3][1][1];
  linEq.matrix[0][3]=xpowsum[0][3];
  linEq.matrix[0][4]=xxpowsum[0][2][1][1];
  linEq.matrix[0][5]=xpowsum[0][2];
  
  linEq.matrix[1][2]=xxpowsum[0][1][1][3];
  linEq.matrix[1][3]=xxpowsum[0][1][1][2];
  linEq.matrix[1][4]=xpowsum[1][3];
  linEq.matrix[1][5]=xpowsum[1][2];
  
  linEq.matrix[2][2]=xxpowsum[0][2][1][2];
  linEq.matrix[2][3]=xxpowsum[0][2][1][1];
  linEq.matrix[2][4]=xxpowsum[0][1][1][2];
  linEq.matrix[2][5]=xxpowsum[0][1][1][1];
  
  linEq.matrix[3][3]=xpowsum[0][2];
  linEq.matrix[3][4]=xxpowsum[0][1][1][1];
  linEq.matrix[3][5]=xpowsum[0][1];
  
  linEq.matrix[4][4]=xpowsum[1][2];
  linEq.matrix[4][5]=xpowsum[1][1];
      
  linEq.matrix[5][5]=xpowsum[0][0];//bottom right entry
  
  //mirror the matrix (top right half mirrored to bottom left half)
  for(i=1;i<linEq.dim;i++)
    for(j=0;j<i;j++)
      linEq.matrix[i][j]=linEq.matrix[j][i];
  
  for(i=0;i<2;i++)
    linEq.vector[i]=mxpowsum[i][2];
  linEq.vector[2]=mxxpowsum[0][1][1][1];
  for(i=3;i<5;i++)
    linEq.vector[i]=mxpowsum[i-3][1];
  linEq.vector[5]=msum;
    
  //solve system of equations and assign values
  if(!(solve_lin_eq(&linEq)==1))
    {
      printf("ERROR: Could not determine fit parameters.\n");
      exit(-1);
    }
  
  //save fit parameters  
  for(i=0;i<linEq.dim;i++)
    a[i]=linEq.solution[i];
    
  //now that the fit is performed, use the fit parameters (and the derivative of the fitting function) to find the minimum
  linEq.dim=2;
  linEq.matrix[0][0]=2*a[0];
  linEq.matrix[0][1]=a[2];
  linEq.matrix[1][0]=a[2];
  linEq.matrix[1][1]=2*a[1];
      
  linEq.vector[0]=-1*a[3];
  linEq.vector[1]=-1*a[4];
  
  //solve system of equations and assign values
  if(!(solve_lin_eq(&linEq)==1))
    {
      printf("ERROR: Could not determine paraboloid center point.\n");
      exit(-1);
    }
  
  //print results
  printf("\nFIT RESULTS\n-----------\n");
  printf("Fit function: f(x,y) = a1*x^2 + a2*y^2 + a3*x*y\n                     + a4*x + a5*y + a6\n\n");
  printf("Coefficients from fit: a1 = %LE\n",a[0]);
  for(i=1;i<6;i++)
    printf("                       a%i = %LE\n",i+1,a[i]);
  printf("\n");
  
  if(a[0]>=0)
    printf("Minimum in x direction, ");
  else
    printf("Maximum in x direction, ");
  printf("x0 = %LE\n",linEq.solution[0]);
  if(a[1]>=0)
    printf("Minimum in y direction, ");
  else
    printf("Maximum in y direction, ");
  printf("y0 = %LE\n",linEq.solution[1]);
  
}
