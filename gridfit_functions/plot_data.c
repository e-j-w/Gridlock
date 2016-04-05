//generate data to be plotted, by selecting datapoints nearest to the fit vertex
void getPlotDataNearMin(const data * d, const parameters * p, const fit_results * fr, plot_data * pd)
{
  int i,j,k;
  long double minDist,mdv;
  for(i=0;i<p->numVar;i++)
    {
      minDist=BIG_NUMBER;
      for(j=0;j<d->lines;j++)
        {
          mdv=(d->x[i][j]) - (fr->fitVert[i]);
          if(mdv<0) mdv=mdv*-1.0; //for some reason abs() doesn't work with long double?
          if(mdv < minDist)
            {
              minDist=mdv;
              pd->fixedParVal[i]=d->x[i][j];
            }
        }
    }
  
  for(i=0;i<p->numVar;i++)
    printf("Parameter %i data point value closest to vertex: %Lf\n",i,pd->fixedParVal[i]);
    
  
  //generate plot data
  int useDataPoint=0;
  if((p->numVar==3)&&(strcmp(p->plotMode,"1d")==0))
    {
      memset(pd->plotDataSize,0,sizeof(pd->plotDataSize));
      for(i=0;i<3;i++)//plot index (x,y,z)
        for(j=0;j<d->lines;j++)
          {
            //check whether the other (non plot index) variables are all at their fixed values
            useDataPoint=1;
            for(k=0;k<3;k++)//parameter index (x,y,z)
              if(k!=i)//is variable fixed?
                if(d->x[k][j]!=pd->fixedParVal[k])
                  useDataPoint=0;
             
            if(useDataPoint==1)
              {
                for(k=0;k<=p->numVar;k++)//parameter index (x,y,z)
                  pd->x[i][k][pd->plotDataSize[i]]=((double)d->x[k][j]);
                pd->plotDataSize[i]++;
              }
          }
      
      for(i=0;i<3;i++)
        printf("Data points in plot %i: %i\n",i,pd->plotDataSize[i]);
    }
  else if((p->numVar==3)&&(strcmp(p->plotMode,"2d")==0))
    {
      memset(pd->plotDataSize,0,sizeof(pd->plotDataSize));
      for(i=0;i<3;i++)//plot index (yz,xz,xy)
        for(j=0;j<d->lines;j++)
          if(d->x[i][j]==pd->fixedParVal[i])
            {
              for(k=0;k<=p->numVar;k++)//parameter index (x,y,z,value)
                pd->x[i][k][pd->plotDataSize[i]]=((double)d->x[k][j]);
              pd->plotDataSize[i]++;
            }
      
      for(i=0;i<3;i++)
        printf("Data points in plot %i: %i\n",i,pd->plotDataSize[i]);
    }
  else if (p->plotData!=0)
    {
      printf("ERROR: Plotting mode not recognized!\nPlot mode: %s\n",p->plotMode);
    }

}

void plotData(const data * d, const parameters * p, const fit_results * fr, plot_data * pd)
{
  int i;
  char str[256];
  plotOpen=1; 
  handle=gnuplot_init();
  printf("\nDATA PLOTS\n----------\nUse 'l' in the plotting window to switch between linear and logarithmic scale.\n");
  
  if((p->numVar==3)&&(strcmp(p->plotMode,"1d")==0))
    {
      for(i=0;i<p->numVar;i++)
        {
          //gnuplot_setstyle(handle,"steps");
          gnuplot_cmd(handle,"set ylabel 'Value'");
          sprintf(str,"set xlabel 'Parameter %i'",i+1);
          gnuplot_cmd(handle,str);
          gnuplot_plot_xy(handle, pd->x[i][i], pd->x[i][p->numVar], pd->plotDataSize[i], "Data");
          printf("Showing plot for parameter %i.  Press [ENTER] to continue.",i);
          getc(stdin);
          gnuplot_resetplot(handle);
        }
    }
  

  return;

}

//function run after CTRL-C, used to clean up temporary files generated
//by the plotting library
void sigint_cleanup()
{
  if(plotOpen==1)
    gnuplot_close(handle); //cleans up temporary files  
  exit(1); 
}
