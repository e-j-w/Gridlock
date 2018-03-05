//filters data assuming it is linearly distributed, by removing data points
//where x/y differs by a value greater than sigma standard deviations
//away from the mean value
void linearFilter(data * d, parameters * p)
{
  int i;
  double avg,stdev;
  
  if(p->numVar>1)
    {
      printf("ERROR: Cannot apply linear filter to data with more than one variable.\n");
      exit(-1);
    }
  
  //get the average difference between x and y
  avg=0.;
  for(i=0;i<d->lines;i++)
    {
      if(d->x[1][i]!=0.)
        avg+=(d->x[0][i]/d->x[1][i]);
    }
  avg/=1.*d->lines;
  
  //get the standard deviation of the difference between x and y
  stdev=0.;
  for(i=0;i<d->lines;i++)
    {
      if(d->x[1][i]!=0.)
        stdev+=((d->x[0][i]/d->x[1][i]) - avg)*((d->x[0][i]/d->x[1][i]) - avg);
    }
  stdev/=(1.0*d->lines - 1.0);
  stdev=sqrt(stdev);
  
  if(p->verbose<1)
    {
      printf("Applying linear filter to data.\n");
      printf("Average: %.3E, Standard Deviation: %.3E, Filter Sigma: %.3E.\n",avg,stdev,p->filterSigma);
    }
  
  //allocate temp data structure
  data *temp=(data*)calloc(1,sizeof(data));
  memcpy(temp,d,sizeof(data));
  
  //filter data
  temp->lines=0;
  for(i=0;i<d->lines;i++)
    if(d->x[1][i]!=0.)
      if((d->x[0][i]/d->x[1][i]) < (avg+(p->filterSigma*stdev)))
        if((d->x[0][i]/d->x[1][i]) > (avg-(p->filterSigma*stdev)))
          {
  	        temp->x[0][temp->lines]=d->x[0][i];
  	        temp->x[1][temp->lines]=d->x[1][i];
  	        temp->lines++;
  	      }
  	      
  if(p->verbose<1)
    printf("%i data points filtered out.\n",d->lines-temp->lines);
  
  //copy back temp data
  memcpy(d,temp,sizeof(data));
  //clean up
  free(temp);
  
}
