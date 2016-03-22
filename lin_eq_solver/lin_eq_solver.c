#include "lin_eq_solver.h"

int solve_lin_eq(lin_eq_type *lin_eq)
{
  lin_eq_type z;
  long double w;
  int i,n;

  n=lin_eq->dim;

  memcpy(z.matrix,lin_eq->matrix,sizeof(lin_eq->matrix));
  w=det(n,&z);

  if(w==0.) 
	  {
	    //printf("Linear equation set has no solutions\n");
	    return 0;
	  }
  for(i=0;i<n;i++)
    {
      memcpy(z.matrix,lin_eq->matrix,sizeof(lin_eq->matrix));
      memcpy(z.matrix[i],lin_eq->vector,sizeof(lin_eq->vector));
      lin_eq->solution[i]=det(n,&z)/w;
    }
  return 1;
}


long double  det(int m, lin_eq_type *lin_eq)
{
  int i,j;
  long double s;

  if(m==1)
    return lin_eq->matrix[0][0];

  if(lin_eq->matrix[m-1][m-1]==0.)
    {
      j=m-1;
      while(lin_eq->matrix[m-1][j]==0 && j>=0) j--;
      if(j<0) 
	      return 0.;
      else
	      for(i=0;i<m;i++)
      	  {
	          s=lin_eq->matrix[i][m-1];
	          lin_eq->matrix[i][m-1]=lin_eq->matrix[i][j];
	          lin_eq->matrix[i][j]=s;
	        }
    }
  for(j=m-2;j>=0;j--)
    for(i=0;i<m;i++)
      lin_eq->matrix[i][j]-=lin_eq->matrix[i][m-1]/lin_eq->matrix[m-1][m-1]*lin_eq->matrix[m-1][j];
  return lin_eq->matrix[m-1][m-1]*det(m-1,lin_eq);
}

