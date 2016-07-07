//generates the sums that will be used when fitting
void generateSums(data * d,const parameters * p)
{
  
  long double powVal=0.;
  long double w=1.;
  int i,j,k,l,m,n,o,q;
  
  //initialize sums (in case this function is called more than once)
	d->msum=0.;
  memset(d->xpowsum,0,POWSIZE*POWSIZE*sizeof(long double));
  memset(d->mxpowsum,0,POWSIZE*POWSIZE*sizeof(long double));
  memset(d->mxxpowsum,0,POWSIZE*POWSIZE*POWSIZE*POWSIZE*sizeof(long double));
  memset(d->xxpowsum,0,POWSIZE*POWSIZE*POWSIZE*POWSIZE*sizeof(long double));
  memset(d->xxxpowsum,0,POWSIZE*POWSIZE*POWSIZE*POWSIZE*POWSIZE*POWSIZE*sizeof(long double));
  
  for(i=0;i<d->lines;i++)//loop over data points
    {
    	w=d->x[p->numVar+1][i]*d->x[p->numVar+1][i];
      d->msum+=d->x[p->numVar][i]/w;
      for(j=0;j<p->numVar;j++)//loop over free parameters
        {
          powVal=1.;
          for(k=0;k<6;k++)//loop over powers
            {
              d->xpowsum[j][k] += powVal/w;
              d->mxpowsum[j][k] += d->x[p->numVar][i]*powVal/w;
              powVal=powVal*d->x[j][i];
            }
            
          for(k=0;k<p->numVar;k++)//loop over free parameters
            for(l=0;l<6;l++)//loop over powers (corresponding to parameter indexed by j)
              for(m=0;m<6;m++)//loop over powers (corresponding to parameter indexed by k)
                {
                  powVal=1.;
                  for(q=0;q<l;q++)
                    powVal=powVal*d->x[j][i];
                  for(q=0;q<m;q++)
                    powVal=powVal*d->x[k][i];
                  d->xxpowsum[j][l][k][m] += powVal/w;
                  d->mxxpowsum[j][l][k][m] += d->x[p->numVar][i]*powVal/w;
                }
          
          for(k=0;k<p->numVar;k++)//loop over free parameters
            for(l=0;l<p->numVar;l++)//loop over free parameters
              for(m=0;m<3;m++)//loop over powers (corresponding to parameter indexed by j)
                for(n=0;n<3;n++)//loop over powers (corresponding to parameter indexed by k)
                  for(o=0;o<3;o++)//loop over powers (corresponding to parameter indexed by l)
                    {
                      powVal=1.;
                      for(q=0;q<m;q++)
                        powVal=powVal*d->x[j][i];
                      for(q=0;q<n;q++)
                        powVal=powVal*d->x[k][i];
                      for(q=0;q<o;q++)
                        powVal=powVal*d->x[l][i];
                      d->xxxpowsum[j][m][k][n][l][o] += powVal/w;         
                    }
        }
    }
}
