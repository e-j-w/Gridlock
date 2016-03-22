#include "gridfit.h"

int main(int argc, char *argv[])
{

  if(argc!=2)
    {
      printf("gridfit filename\n");
      printf("\nPerforms grid minimization on the data in the specified file.  The file should be in plaintext, with the first 3 columns corresponding to free parameters and the 4th column corresponding to the grid point value.\n");
      exit(-1);
    }

  //int numPar=atoi(argv[2]);
  int numPar=3; //eventually plan to extend to n=2,n=1 case
 
  if((inp=fopen(argv[1],"r"))==NULL)
    {
      printf("\nERROR: input file %s can not be opened.\n",argv[1]);
      exit(-1);
    }

  //initialize values
  //int ndf=0.;
  int lines=0;
  int linenum=0;
  memset(x,0,sizeof(x));
  msum=0.;
  memset(xpowsum,0,sizeof(xpowsum));
  memset(mxpowsum,0,sizeof(mxpowsum));
  memset(xxpowsum,0,sizeof(xxpowsum));

  while(!(feof(inp)))//go until the end of file is reached
    {
      if(fgets(str,256,inp)!=NULL)
        {
          if(sscanf(str,"%Lf %Lf %Lf %Lf %Lf %Lf",&x[0][lines],&x[1][lines],&x[2][lines],&x[3][lines],&x[4][lines],&x[5][lines])==numPar+1)
            {
              lines++;
            }
          else
            {
              printf("WARNING: Improperly formatted data on line %i of input file.\n",linenum+1);
            }
          linenum++;
        }
    }

  fclose(inp);
  printf("%i lines of data read from file: %s\n",lines,argv[1]);


  //construct sums
  long double powVal; 
  for(i=0;i<lines;i++)//loop over data points
    {
      msum+=x[numPar][i];
      for(j=0;j<numPar;j++)//loop over free parameters
        {
          powVal=1.;
          for(k=0;k<5;k++)//loop over powers
            {
              xpowsum[j][k] += powVal;
              mxpowsum[j][k] += x[numPar][i]*powVal;
              powVal=powVal*x[j][i];
            }
            
          for(k=0;k<numPar;k++)//loop over free parameters
            for(l=0;l<5;l++)//loop over powers (corresponding to parameter indexed by j)
              for(m=0;m<5;m++)//loop over powers (corresponding to parameter indexed by k)
                {
                  powVal=1.;
                  for(p=0;p<l;p++)
                    powVal=powVal*x[j][i];
                  for(p=0;p<m;p++)
                    powVal=powVal*x[k][i];
                  xxpowsum[j][l][k][m] += powVal;
                  mxxpowsum[j][l][k][m] += x[numPar][i]*powVal;
                }
          
          for(k=0;k<numPar;k++)//loop over free parameters
            for(l=0;l<numPar;l++)//loop over free parameters
              for(m=0;m<3;m++)//loop over powers (corresponding to parameter indexed by j)
                for(n=0;n<3;n++)//loop over powers (corresponding to parameter indexed by k)
                  for(o=0;o<3;o++)//loop over powers (corresponding to parameter indexed by l)
                    {
                      powVal=1.;
                      for(p=0;p<m;p++)
                        powVal=powVal*x[j][i];
                      for(p=0;p<n;p++)
                        powVal=powVal*x[k][i];
                      for(p=0;p<o;p++)
                        powVal=powVal*x[l][i];
                      xxxpowsum[j][m][k][n][l][o] += powVal;         
                    }
        }
    }
  
  
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
      printf("ERROR: Could not find a proper minimization of chisq!\n");
      exit(-1);
    }
    
  for(i=0;i<linEq.dim;i++)
    a[i]=linEq.solution[i];
  
  printf("\nFIT RESULTS\n-----------\n");
  printf("Coefficients from paraboloid fit:\n");
  for(i=0;i<linEq.dim;i++)
    printf("a%i = %LE\n",i+1,a[i]);
    
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
      printf("ERROR: Could not find a proper minimization of chisq!\n");
      exit(-1);
    }
    
  printf("Paraboloid minimum:\n");
  printf("x0 = %LE\n",linEq.solution[0]);
  printf("y0 = %LE\n",linEq.solution[1]);
  printf("z0 = %LE\n",linEq.solution[2]);

  return 0; //great success
}
