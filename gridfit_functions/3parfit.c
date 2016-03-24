//fit data to a paraboloid of the form
//f(x,y,z) = a1*x^2 + a2*y^2 + a3*z^2 + a4*x*y + a5*x*z + a6*y*z +a7*x + a8*y + a9*z + a10
void fit3Par()
{
  //construct equations (n=3 specific case)
  linEq.dim=10;
  
  for(i=0;i<3;i++)//loop over free parameters
    for(j=0;j<3;j++)//loop over free parameters
      linEq.matrix[i][j]=xxpowsum[i][2][j][2];//top-left 3x3 entries
      
  linEq.matrix[0][3]=xxpowsum[0][3][1][1];
  linEq.matrix[0][4]=xxpowsum[0][3][2][1];
  linEq.matrix[0][5]=xxxpowsum[0][2][1][1][2][1];
  linEq.matrix[0][6]=xpowsum[0][3];
  linEq.matrix[0][7]=xxpowsum[0][2][1][1];
  linEq.matrix[0][8]=xxpowsum[0][2][2][1];
  linEq.matrix[0][9]=xpowsum[0][2];
  
  linEq.matrix[1][3]=xxpowsum[0][1][1][3];
  linEq.matrix[1][4]=xxxpowsum[0][1][1][2][2][1];
  linEq.matrix[1][5]=xxpowsum[1][3][2][1];
  linEq.matrix[1][6]=xxpowsum[0][1][1][2];
  linEq.matrix[1][7]=xpowsum[1][3];
  linEq.matrix[1][8]=xxpowsum[1][2][2][1];
  linEq.matrix[1][9]=xpowsum[1][2];
  
  linEq.matrix[2][3]=xxxpowsum[0][1][1][1][2][2];
  linEq.matrix[2][4]=xxpowsum[0][1][2][3];
  linEq.matrix[2][5]=xxpowsum[1][1][2][3];
  linEq.matrix[2][6]=xxpowsum[0][1][2][2];
  linEq.matrix[2][7]=xxpowsum[1][1][2][2];
  linEq.matrix[2][8]=xpowsum[2][3];
  linEq.matrix[2][9]=xpowsum[2][2];

  linEq.matrix[3][3]=xxpowsum[0][2][1][2];  
  linEq.matrix[3][4]=xxxpowsum[0][2][1][1][2][1];
  linEq.matrix[3][5]=xxxpowsum[0][1][1][2][2][1];
  linEq.matrix[3][6]=xxpowsum[0][2][1][1];
  linEq.matrix[3][7]=xxpowsum[0][1][1][2];
  linEq.matrix[3][8]=xxxpowsum[0][1][1][1][2][1];
  linEq.matrix[3][9]=xxpowsum[0][1][1][1];
  
  linEq.matrix[4][4]=xxpowsum[0][2][2][2];
  linEq.matrix[4][5]=xxxpowsum[0][1][1][1][2][2];
  linEq.matrix[4][6]=xxpowsum[0][2][2][1];
  linEq.matrix[4][7]=xxxpowsum[0][1][1][1][2][1];
  linEq.matrix[4][8]=xxpowsum[0][1][2][2];
  linEq.matrix[4][9]=xxpowsum[0][1][2][1];
  
  linEq.matrix[5][5]=xxpowsum[1][2][2][2];
  linEq.matrix[5][6]=xxxpowsum[0][1][1][1][2][1];
  linEq.matrix[5][7]=xxpowsum[1][2][2][1];
  linEq.matrix[5][8]=xxpowsum[1][1][2][2];
  linEq.matrix[5][9]=xxpowsum[1][1][2][1];
  
  linEq.matrix[6][6]=xpowsum[0][2];
  linEq.matrix[6][7]=xxpowsum[0][1][1][1];
  linEq.matrix[6][8]=xxpowsum[0][1][2][1];
  linEq.matrix[6][9]=xpowsum[0][1];
  
  linEq.matrix[7][7]=xpowsum[1][2];
  linEq.matrix[7][8]=xxpowsum[1][1][2][1];
  linEq.matrix[7][9]=xpowsum[1][1];
  
  linEq.matrix[8][8]=xpowsum[2][2];
  linEq.matrix[8][9]=xpowsum[2][1];
      
  linEq.matrix[9][9]=xpowsum[0][0];//bottom right entry
  
  //mirror the matrix (top right half mirrored to bottom left half)
  for(i=1;i<linEq.dim;i++)
    for(j=0;j<i;j++)
      linEq.matrix[i][j]=linEq.matrix[j][i];
  
  for(i=0;i<3;i++)
    linEq.vector[i]=mxpowsum[i][2];
  linEq.vector[3]=mxxpowsum[0][1][1][1];
  linEq.vector[4]=mxxpowsum[0][1][2][1];
  linEq.vector[5]=mxxpowsum[1][1][2][1];
  for(i=6;i<9;i++)
    linEq.vector[i]=mxpowsum[i-6][1];
  linEq.vector[9]=msum;
    
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
  linEq.dim=3;
  linEq.matrix[0][0]=2*a[0];
  linEq.matrix[0][1]=a[3];
  linEq.matrix[0][2]=a[4];
  linEq.matrix[1][1]=2*a[1];
  linEq.matrix[1][2]=a[5];
  linEq.matrix[2][2]=2*a[2];
  //mirror the matrix (top right half mirrored to bottom left half)
  for(i=1;i<linEq.dim;i++)
    for(j=0;j<i;j++)
      linEq.matrix[i][j]=linEq.matrix[j][i];
      
  linEq.vector[0]=-1*a[6];
  linEq.vector[1]=-1*a[7];
  linEq.vector[2]=-1*a[8];
  
  //solve system of equations and assign values
  if(!(solve_lin_eq(&linEq)==1))
    {
      printf("ERROR: Could not determine paraboloid center point.\n");
      exit(-1);
    }
  
  //print results
  printf("\nFIT RESULTS\n-----------\n");
  printf("Fit function: f(x,y,z) = a1*x^2 + a2*y^2 + a3*z^2\n                       + a4*x*y + a5*x*z + a6*y*z\n                       + a7*x + a8*y + a9*z + a10\n\n");
  printf("Coefficients from fit: a1 = %LE\n",a[0]);
  for(i=1;i<10;i++)
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
  if(a[2]>=0)
    printf("Minimum in z direction, ");
  else
    printf("Maximum in z direction, ");  
  printf("z0 = %LE\n",linEq.solution[2]);
  
  long double fitVal=a[0]*linEq.solution[0]*linEq.solution[0] + a[1]*linEq.solution[1]*linEq.solution[1] + a[2]*linEq.solution[2]*linEq.solution[2] + a[3]*linEq.solution[0]*linEq.solution[1] + a[4]*linEq.solution[0]*linEq.solution[2] + a[5]*linEq.solution[1]*linEq.solution[2] + a[6]*linEq.solution[0] + a[7]*linEq.solution[1] + a[8]*linEq.solution[2] + a[9];
  
  printf("\nf(x0,y0,z0) = %LE\n",fitVal); 
  
}
