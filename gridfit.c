//definitions
#include "gridfit.h"
//functions
#include "2parfit.c"
#include "3parfit.c"

int main(int argc, char *argv[])
{

  if(argc!=3)
    {
      printf("\ngridfit filename n\n");
      printf("\nPerforms grid minimization on the data in the specified file.  The file should be in plaintext, with the first n columns corresponding to free parameters and the (n+1)th column corresponding to the grid point value.\n\n");
      exit(-1);
    }

  int numPar=atoi(argv[2]);
  if((numPar<2)||(numPar>3))
    {
      printf("ERROR: the number of free parameters n must be 2 or 3.\nAborting...\n");
      exit(-1);
    }
 
  if((inp=fopen(argv[1],"r"))==NULL)
    {
      printf("\nERROR: input file %s can not be opened.\n",argv[1]);
      exit(-1);
    }

  //initialize values
  int lines=0;
  int linenum=0;
  memset(x,0,sizeof(x));
  msum=0.;
  memset(xpowsum,0,sizeof(xpowsum));
  memset(mxpowsum,0,sizeof(mxpowsum));
  memset(mxxpowsum,0,sizeof(mxxpowsum));
  memset(xxpowsum,0,sizeof(xxpowsum));
  memset(xxxpowsum,0,sizeof(xxxpowsum));

  //import data from file
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
  
  
  //call specific fitting routines depending on the number of free parameters
  if(numPar==2)
    fit2Par();
  else if(numPar==3)
    fit3Par();
    
  return 0; //great success
}
