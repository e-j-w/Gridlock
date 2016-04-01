//prints info and statistics for data
void printDataInfo(const data * d, const parameters * p)
{
  int i,j;
  int maxInd[NUM_LIST],minInd[NUM_LIST];
  int numMax=0;
  int numMin=0;
  long double maxVal=0.;
  long double minVal=BIG_NUMBER;
  for(i=0;i<d->lines;i++)
    {
      if(d->x[p->numVar][i]>maxVal)
        {
          maxVal=d->x[p->numVar][i];
          for(j=1;j<NUM_LIST;j++)
            maxInd[NUM_LIST-j]=maxInd[NUM_LIST-(j+1)];//shift values down the array
          maxInd[0]=i;//record the new max index
          numMax++;
        }
      if(d->x[p->numVar][i]<minVal)
        {
          minVal=d->x[p->numVar][i];
          for(j=1;j<NUM_LIST;j++)
            minInd[NUM_LIST-j]=minInd[NUM_LIST-(j+1)];//shift values down the array
          minInd[0]=i;//record the new min index
          numMin++;
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
