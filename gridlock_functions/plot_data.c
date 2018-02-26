//forward declarations
long double evalLin(long double, const fit_results *);
long double evalPoly3(long double, const fit_results *);
long double eval2ParPoly3(long double,long double, const fit_results *);
long double eval1Par(long double, const fit_results *);
long double eval2Par(long double,long double, const fit_results *);
long double eval3Par(long double,long double,long double, const fit_results *);

//generate data to be plotted
//for multidimensional paraboloid fits, do this by by selecting datapoints nearest to the fit vertex
void preparePlotData(const data * d, const parameters * p, const fit_results * fr, plot_data * pd)
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
  
  //generate plot data
  int useDataPoint=0;
  if(strcmp(p->plotMode,"1d")==0)
    {
      pd->numPlots=p->numVar;
      memset(pd->plotDataSize,0,sizeof(pd->plotDataSize));
      for(i=0;i<pd->numPlots;i++)//plot index (x,y,z)
        for(j=0;j<d->lines;j++)
          {
            //check whether the other (non plot index) variables are all at their fixed values
            useDataPoint=1;
            for(k=0;k<p->numVar;k++)//parameter index (x,y,z)
              if(k!=i)//is variable fixed?
                if(d->x[k][j]!=pd->fixedParVal[k])
                  useDataPoint=0;
             
            if(useDataPoint==1)
              {
                for(k=0;k<=p->numVar;k++)//parameter index (x,y,z)
                  pd->data[i][k][pd->plotDataSize[i]]=((double)d->x[k][j]);
                  
                pd->plotDataSize[i]++;
              }
          }
    }
  else if((p->numVar==3)&&(strcmp(p->plotMode,"2d")==0))
    {
      pd->numPlots=3;
      memset(pd->plotDataSize,0,sizeof(pd->plotDataSize));
      for(i=0;i<pd->numPlots;i++)//plot index (yz,xz,xy)
        for(j=0;j<d->lines;j++)
          if(d->x[i][j]==pd->fixedParVal[i])
            {
              for(k=0;k<=p->numVar;k++)//parameter index (x,y,z,value)
                pd->data[i][k][pd->plotDataSize[i]]=((double)d->x[k][j]);

              pd->plotDataSize[i]++;
            }
    }
  else if((p->numVar==2)&&(strcmp(p->plotMode,"2d")==0))
    {
      pd->numPlots=1;
      memset(pd->plotDataSize,0,sizeof(pd->plotDataSize));
      for(i=0;i<d->lines;i++)
        {
          //copy over data to plot
          for(j=0;j<=p->numVar;j++)//parameter index (x,y,z,value)
            pd->data[0][j][pd->plotDataSize[0]]=((double)d->x[j][i]);
            
          pd->plotDataSize[0]++;
        }
    }
  else if((p->numVar==3)&&(strcmp(p->plotMode,"3d")==0))
    {
      pd->numPlots=1;
      memset(pd->plotDataSize,0,sizeof(pd->plotDataSize));
      for(i=0;i<d->lines;i++)
        {
          //copy over data to plot
          for(j=0;j<=p->numVar;j++)//parameter index (x,y,z,value)
            pd->data[0][j][pd->plotDataSize[0]]=((double)d->x[j][i]);

          pd->plotDataSize[0]++;
        }
    }
  else if (p->plotData!=0)
    {
      printf("ERROR: Plotting mode '%s' is not availiable for the fit type used (%s).\n",p->plotMode,p->fitType);
      if(p->numVar==1)
        printf("Available plot modes: 1d.\n");
      else if(p->numVar==2)
        printf("Available plot modes: 1d, 2d.\n");
      else if(p->numVar==3)
        printf("Available plot modes: 1d, 2d, 3d.\n");
      else
        printf("No plot modes available for this fit type.\n");
      exit(-1);
    }

  //determine whether or not to use scientific notation for labels
  for(i=0;i<pd->numPlots;i++)
    for(j=0;j<p->numVar;j++)
      for(k=0;k<pd->plotDataSize[i];k++)
        if((fabs(pd->data[i][j][k])<0.001)&&(pd->data[i][j][k]!=0.))
          pd->axisLabelStyle[i][j]=1;
        else
        	{
        		pd->axisLabelStyle[i][j]=0;
        		break;//if any data point is not small, use regular labels
        	}
  
  //copy over minimum/maximum values (for determining plotting ranges)
  pd->min_m=(double)d->min_m;
  pd->max_m=(double)d->max_m;

  //generate fit plot data
  pd->numFitPtsPerVar=20;
  pd->numFitPlotPts=1;
  double floorFactor;
  for(i=0;i<p->numVar;i++)
    pd->numFitPlotPts=pd->numFitPlotPts*pd->numFitPtsPerVar; //compute numFitPtsPerVar^numVar

  for(i=0;i<pd->numPlots;i++)
    for(j=0;j<pd->numFitPlotPts;j++)
      if(j<MAXFILELENGTH)
        {
          floorFactor=1;
          for(k=0;k<p->numVar;k++)
            {
              if(k!=0)
                floorFactor=floorFactor*pd->numFitPtsPerVar;
              
              //variable not fixed
              //had to use a spreadsheet to visualize this haHAA
              pd->fit[i][k][j]=d->min_x[k] + (((int)(floor((double)j/floorFactor))%pd->numFitPtsPerVar)/(double)pd->numFitPtsPerVar)*(d->max_x[k] - d->min_x[k]);

              //handle fixed variable cases
              if(strcmp(p->plotMode,"1d")==0)
                {
                  if(k!=i)//variable fixed
                    pd->fit[i][k][j]=(double)pd->fixedParVal[k];
                }
              else if((p->numVar==3)&&(strcmp(p->plotMode,"2d")==0))
                {
                  if(k==i)//variable fixed
                    pd->fit[i][k][j]=(double)pd->fixedParVal[k];
                }
            }
          if(strcmp(p->fitType,"par1")==0)
            pd->fit[i][p->numVar][j]=(double)eval1Par(pd->fit[i][0][j],fr);
          else if((strcmp(p->fitType,"lin")==0)||(strcmp(p->fitType,"lin_deming")==0))
            pd->fit[i][p->numVar][j]=(double)evalLin(pd->fit[i][0][j],fr);
          else if(strcmp(p->fitType,"poly3")==0)
            pd->fit[i][p->numVar][j]=(double)evalPoly3(pd->fit[i][0][j],fr);
          else if(strcmp(p->fitType,"par2")==0)
            pd->fit[i][p->numVar][j]=((double)eval2Par(pd->fit[i][0][j],pd->fit[i][1][j],fr));
          else if(strcmp(p->fitType,"2parpoly3")==0)
            pd->fit[i][p->numVar][j]=((double)eval2ParPoly3(pd->fit[i][0][j],pd->fit[i][1][j],fr));
          else if(strcmp(p->fitType,"par3")==0)
            pd->fit[i][p->numVar][j]=(double)eval3Par(pd->fit[i][0][j],pd->fit[i][1][j],pd->fit[i][2][j],fr);
          else
            printf("WARNING: Unknown fit type '%s', cannot plot fit.\n",p->fitType);
        }
    

}

//handles the gnuplot prompt
void plotPrompt(int cont)
{
	int c;
	char inp[256];
	if(cont==0)
 		printf("Enter 'g' for a gnuplot prompt or press [ENTER] to exit. ");
	else
		printf("Enter 'g' for a gnuplot prompt or press [ENTER] to continue. ");
	c=getc(stdin);
	if(c=='g')
		{
			printf("Enter 'exit' to return from the gnuplot prompt.\n");
			fgets(inp,256,stdin);
			while(strcmp(inp,"exit\n")!=0)
				{
					gnuplot_cmd(handle,inp);
					printf("gnuplot > ");
					fgets(inp,256,stdin);
				}
		}
	return;
}

void plotData(const parameters * p, fit_results * fr, plot_data * pd)
{
  int i;
  char * str=(char*)calloc(256,sizeof(char));
  plotOpen=1; 
  handle=gnuplot_init();
    
  printf("\nDATA PLOTS\n----------\n");
  
  if(strcmp(p->plotMode,"1d")==0)
    {
      for(i=0;i<p->numVar;i++)
        {
          gnuplot_setstyle(handle,"points"); //set style for grid points
          if(strcmp(p->dataType,"chisq")==0)
            gnuplot_cmd(handle,"set ylabel 'Chisq'");
          else  
            gnuplot_cmd(handle,"set ylabel 'Value'");
          sprintf(str,"set xlabel 'Parameter %i'",i+1);
          gnuplot_cmd(handle,str);
          gnuplot_plot_xy(handle, pd->data[i][i], pd->data[i][p->numVar], pd->plotDataSize[i], "Data");
          if(pd->axisLabelStyle[i][i]==1)
            gnuplot_cmd(handle,"set format x '%%12.2E'");
          if(pd->axisLabelStyle[i][p->numVar]==1)
            gnuplot_cmd(handle,"set format y '%%12.2E'");
          gnuplot_setstyle(handle,"lines");//set style for fit data
          if(p->numVar==1)
            {
              gnuplot_plot_xy(handle, pd->fit[i][i], pd->fit[i][p->numVar], pd->numFitPlotPts, "Fit");
            }  
          else
            {
              if(i==0)
                gnuplot_plot_xy(handle, pd->fit[i][i], pd->fit[i][p->numVar], pd->numFitPtsPerVar, "Fit");
              else if(i==1)
                gnuplot_plot_xygrid(handle, pd->fit[i][i], pd->fit[i][p->numVar], pd->numFitPtsPerVar*pd->numFitPtsPerVar, pd->numFitPtsPerVar, pd->numFitPtsPerVar, 0, "Fit");
              else if(i==2)
                gnuplot_plot_xygrid(handle, pd->fit[i][i], pd->fit[i][p->numVar], pd->numFitPlotPts, pd->numFitPtsPerVar, pd->numFitPtsPerVar*pd->numFitPtsPerVar, 0, "Fit");
            }
          //strcpy(str,fr->fitForm[i]);//retrieve fit data functional form
          //gnuplot_plot_equation(handle, str, "Fit (function)");
          //plot confidence intervals
          if(p->plotCI==1)
          	{
          		gnuplot_plot_xy(handle,fr->ciXVal[i],fr->ciUVal[i],CI_DIM,"Upper 1-sigma confidence band");
          		gnuplot_plot_xy(handle,fr->ciXVal[i],fr->ciLVal[i],CI_DIM,"Lower 1-sigma confidence band");
          		/*strcpy(str,fr->ciUForm[0]);
          		gnuplot_plot_equation(handle, str, "Upper 95% confidence band");
          		strcpy(str,fr->ciLForm[0]);
          		gnuplot_plot_equation(handle, str, "Lower 95% confidence band");
          		strcpy(str,fr->piUForm[0]);
          		gnuplot_plot_equation(handle, str, "Upper 95% prediction band");
          		strcpy(str,fr->piLForm[0]);
          		gnuplot_plot_equation(handle, str, "Lower 95% prediction band");*/
          	}
          printf("Showing plot for parameter %i.\n",i+1);
          if(p->numVar==3)
            {
              if(i==0)
                printf("Parameter %i fixed to %Lf\nParameter %i fixed to %Lf\n",2,pd->fixedParVal[1],3,pd->fixedParVal[2]);
              if(i==1)
                printf("Parameter %i fixed to %Lf\nParameter %i fixed to %Lf\n",1,pd->fixedParVal[0],3,pd->fixedParVal[2]);
              if(i==2)
                printf("Parameter %i fixed to %Lf\nParameter %i fixed to %Lf\n",1,pd->fixedParVal[0],2,pd->fixedParVal[1]);
            }
          else if(p->numVar==2)
            {
              if(i==0)
                printf("Parameter %i fixed to %Lf\n",2,pd->fixedParVal[1]);
              if(i==1)
                printf("Parameter %i fixed to %Lf\n",1,pd->fixedParVal[0]);
            }
          printf("%i data points available for plot.\n",pd->plotDataSize[i]);
          if(i<(p->numVar-1))
            plotPrompt(1);
          else
            plotPrompt(0);
          gnuplot_resetplot(handle);
        }
    }
  else if(strcmp(p->plotMode,"2d")==0)
    {
      if(p->numVar==3)
        {
          for(i=0;i<p->numVar;i++)
            {
              gnuplot_setstyle(handle,"points"); //set style for grid points
              if(i==0)
                {
                  gnuplot_plot_xyz(handle, pd->data[i][1], pd->data[i][2], pd->data[i][p->numVar], pd->plotDataSize[i], "Data");
                  if(pd->axisLabelStyle[i][1]==1)
                    gnuplot_cmd(handle,"set format x '%%12.2E'");
                  if(pd->axisLabelStyle[i][2]==1)
                    gnuplot_cmd(handle,"set format y '%%12.2E'");
                  if(pd->axisLabelStyle[i][p->numVar]==1)
                    gnuplot_cmd(handle,"set format z '%%12.2E'");
                  sprintf(str,"set xlabel 'Parameter 2'; set ylabel 'Parameter 3'");
                }
              else if(i==1)
                {
                  gnuplot_plot_xyz(handle, pd->data[i][0], pd->data[i][2], pd->data[i][p->numVar], pd->plotDataSize[i], "Data");
                  if(pd->axisLabelStyle[i][0]==1)
                    gnuplot_cmd(handle,"set format x '%%12.2E'");
                  if(pd->axisLabelStyle[i][2]==1)
                    gnuplot_cmd(handle,"set format y '%%12.2E'");
                  if(pd->axisLabelStyle[i][p->numVar]==1)
                    gnuplot_cmd(handle,"set format z '%%12.2E'");
                  sprintf(str,"set xlabel 'Parameter 1'; set ylabel 'Parameter 3'");
                }
              else if(i==2)
                {
                  gnuplot_plot_xyz(handle, pd->data[i][0], pd->data[i][1], pd->data[i][p->numVar], pd->plotDataSize[i], "Data");
                  if(pd->axisLabelStyle[i][0]==1)
                    gnuplot_cmd(handle,"set format x '%%12.2E'");
                  if(pd->axisLabelStyle[i][1]==1)
                    gnuplot_cmd(handle,"set format y '%%12.2E'");
                  if(pd->axisLabelStyle[i][p->numVar]==1)
                    gnuplot_cmd(handle,"set format z '%%12.2E'");
                  sprintf(str,"set xlabel 'Parameter 1'; set ylabel 'Parameter 2'");
                }
              if(strcmp(p->dataType,"chisq")==0)
                gnuplot_cmd(handle,"set zlabel 'Chisq'");
              else  
                gnuplot_cmd(handle,"set zlabel 'Value'");
              gnuplot_cmd(handle,str);
              gnuplot_setstyle(handle,"lines");
              gnuplot_cmd(handle,"set grid");//set style for fit data
              if(i==0)
                gnuplot_plot_xyzgrid(handle, pd->fit[i][1], pd->fit[i][2], pd->fit[i][p->numVar], pd->numFitPlotPts, pd->numFitPtsPerVar*pd->numFitPtsPerVar, pd->numFitPtsPerVar, 0, "Fit");
              else if(i==1)
                gnuplot_plot_xyzgrid(handle, pd->fit[i][0], pd->fit[i][2], pd->fit[i][p->numVar], pd->numFitPlotPts, pd->numFitPtsPerVar, 0, pd->numFitPtsPerVar, "Fit");
              else if(i==2)
                gnuplot_plot_xyzgrid(handle, pd->fit[i][0], pd->fit[i][1], pd->fit[i][p->numVar], pd->numFitPtsPerVar*pd->numFitPtsPerVar, pd->numFitPtsPerVar, 0, 0, "Fit");
              //strcpy(str,fr->fitForm[i]);//retrieve fit data functional form
              //gnuplot_plot_equation(handle, str, "Fit (function)");
              printf("Showing surface plot with parameter %i fixed to %Lf\n",i+1,pd->fixedParVal[i]);
              printf("%i data points available for plot.\n",pd->plotDataSize[i]);
              if(i<(p->numVar-1))
                plotPrompt(1);
              else
                plotPrompt(0);
              gnuplot_resetplot(handle);
            }
        }
      else if(p->numVar==2)
        {
          gnuplot_setstyle(handle,"points"); //set style for grid points
          gnuplot_plot_xyz(handle, pd->data[0][0], pd->data[0][1], pd->data[0][p->numVar], pd->plotDataSize[0], "Data");
          if(pd->axisLabelStyle[0][0]==1)
            gnuplot_cmd(handle,"set format x '%%12.2E'");
          if(pd->axisLabelStyle[0][1]==1)
            gnuplot_cmd(handle,"set format y '%%12.2E'");
          if(pd->axisLabelStyle[0][p->numVar]==1)
            gnuplot_cmd(handle,"set format z '%%12.2E'"); 
          if(strcmp(p->dataType,"chisq")==0)
            gnuplot_cmd(handle,"set zlabel 'Chisq'");
          else  
            gnuplot_cmd(handle,"set zlabel 'Value'");
          sprintf(str,"set xlabel 'Parameter 1'; set ylabel 'Parameter 2'");
          gnuplot_cmd(handle,str);
          gnuplot_setstyle(handle,"lines");
          gnuplot_cmd(handle,"set grid");//set style for fit data
          //gnuplot_cmd(handle,"set dgrid3d 30,30 qnorm 2");//set style for fit data
          gnuplot_plot_xyzgrid(handle, pd->fit[0][0], pd->fit[0][1], pd->fit[0][p->numVar], pd->numFitPlotPts, pd->numFitPtsPerVar, 0, 0, "Fit");
          //strcpy(str,fr->fitForm[0]);//retrieve fit data functional form
          //gnuplot_plot_equation(handle, str, "Fit (function)");
          printf("Showing surface plot.\n");
          printf("%i data points available for plot.\n",pd->plotDataSize[0]);
          plotPrompt(0);
          gnuplot_resetplot(handle);
        }
    }
  else if(strcmp(p->plotMode,"3d")==0)
    {
      if(p->numVar==3)
        {
          gnuplot_setstyle(handle,"points"); //set style for grid points
          gnuplot_plot_xyza(handle, pd->data[0][0], pd->data[0][1], pd->data[0][2], pd->data[0][p->numVar], pd->plotDataSize[0], "Data");
          if(pd->axisLabelStyle[0][0]==1)
            gnuplot_cmd(handle,"set format x '%%12.2E'");
          if(pd->axisLabelStyle[0][1]==1)
            gnuplot_cmd(handle,"set format y '%%12.2E'");
          if(pd->axisLabelStyle[0][2]==1)
            gnuplot_cmd(handle,"set format z '%%12.2E'");
          sprintf(str,"set xlabel 'Parameter 1'; set ylabel 'Parameter 2'; set zlabel 'Parameter 3'");
          gnuplot_cmd(handle,str);
          sprintf(str,"set cbrange [%f:%f]",pd->min_m,pd->max_m); //set the color bar range
          gnuplot_cmd(handle,str);
          if(strcmp(p->fitType,"par3")==0)
            {
              gnuplot_setcolor(handle,"black");
              gnuplot_cmd(handle,"set pointsize 1.5");
              double xVert=(double)fr->fitVert[0];
              double yVert=(double)fr->fitVert[1];
              double zVert=(double)fr->fitVert[2];
              gnuplot_plot_xyz(handle, &xVert, &yVert, &zVert, 1, "Fit Vertex");
            }
          printf("Showing heatmap plot.\n");
          printf("%i data points available for plot.\n",pd->plotDataSize[0]);
          plotPrompt(0);
          gnuplot_resetplot(handle);
        }
    }
  else
    {
      printf("ERROR: The plot mode '%s' defined in the data file '%s' is not supported!\n",p->plotMode,p->filename);
      exit(-1);
    }

}

//function run after CTRL-C, used to clean up temporary files generated
//by the plotting library
void sigint_cleanup()
{
  if(plotOpen==1)
    gnuplot_close(handle); //cleans up temporary files  
  exit(1); 
}
