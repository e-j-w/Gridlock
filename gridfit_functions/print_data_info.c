//prints info and statistics for data
void printDataInfo(const data * d, const parameters * p)
{
  int i,j,k;
  int maxInd[NUM_LIST],minInd[NUM_LIST];
  long double maxVal[NUM_LIST],minVal[NUM_LIST];
  for(i=0;i<NUM_LIST;i++)
    {
      maxVal[i]=0.;
      minVal[i]=BIG_NUMBER;
    }
  int numMax=0;
  int numMin=0;
  for(i=0;i<d->lines;i++)
    {
      for(j=0;j<NUM_LIST;j++)
        if(d->x[p->numVar][i]>maxVal[j])
          { 
            for(k=1;k<NUM_LIST-j;k++)
              {
                maxInd[NUM_LIST-k]=maxInd[NUM_LIST-(k+1)];//shift values down the array
                maxVal[NUM_LIST-k]=maxVal[NUM_LIST-(k+1)];//shift values down the array
              }
            maxVal[j]=d->x[p->numVar][i];
            maxInd[j]=i;//record the new max index entry
            numMax++;
            break;
          }
      for(j=0;j<NUM_LIST;j++)
        if(d->x[p->numVar][i]<minVal[j])
          {
            for(k=1;k<NUM_LIST-j;k++)
              {
                minInd[NUM_LIST-k]=minInd[NUM_LIST-(k+1)];//shift values down the array
                minVal[NUM_LIST-k]=minVal[NUM_LIST-(k+1)];//shift values down the array
              }
            minVal[j]=d->x[p->numVar][i];
            minInd[j]=i;//record the new min index entry
            numMin++;
            break;
          }
    }
  printf("\nData minimum value(s): %0.3LE at [",d->x[p->numVar][minInd[0]]);
  for(i=0;i<p->numVar;i++)
    printf(" %0.3LE ",d->x[i][minInd[0]]);
  printf("]\n");
  for(i=1;i<numMin;i++)
    if(i<NUM_LIST)
      {
        printf("                       %0.3LE at [",d->x[p->numVar][minInd[i]]);
        for(j=0;j<p->numVar;j++)
          printf(" %0.3LE ",d->x[j][minInd[i]]);
        printf("]\n");
      }

  printf("\nData maximum value(s): %0.3LE at [",d->x[p->numVar][maxInd[0]]);
  for(i=0;i<p->numVar;i++)
    printf(" %0.3LE ",d->x[i][maxInd[0]]);
  printf("]\n");
  for(i=1;i<numMax;i++)
    if(i<NUM_LIST)
      {
        printf("                       %0.3LE at [",d->x[p->numVar][maxInd[i]]);
        for(j=0;j<p->numVar;j++)
          printf(" %0.3LE ",d->x[j][maxInd[i]]);
        printf("]\n");
      }
  
}
