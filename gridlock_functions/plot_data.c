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
        if((abs(pd->data[i][j][k])<0.001)&&(pd->data[i][j][k]!=0.))
          pd->axisLabelStyle[i][j]=1;
        else
        	{
        		pd->axisLabelStyle[i][j]=0;
        		break;//if any data point is not small, use regular labels
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
          //generate fit data functional forms
          if(strcmp(p->fitType,"par3")==0)
            str=plotForm3Par(p,fr,pd,i);
          else if(strcmp(p->fitType,"par2")==0)
            str=plotForm2Par(p,fr,pd,i);
          else if(strcmp(p->fitType,"par1")==0)
            str=plotForm1Par(p,fr,pd,i);
          else if((strcmp(p->fitType,"lin")==0)||(strcmp(p->fitType,"lin_deming")==0))
            strcpy(str,fr->fitForm[0]);
          else if(strcmp(p->fitType,"poly3")==0)
            str=plotFormPoly3(p,fr,pd,i);
          else if(strcmp(p->fitType,"2parpoly3")==0)
            str=plotForm2ParPoly3(p,fr,pd,i);
          gnuplot_plot_equation(handle, str, "Fit");
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
              if(i==1)
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
              if(i==2)
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
              //generate fit data functional forms
              if(strcmp(p->fitType,"par3")==0)
                str=plotForm3Par(p,fr,pd,i);
              gnuplot_plot_equation(handle, str, "Fit");
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
          //generate fit data functional forms
          if(strcmp(p->fitType,"par2")==0)
            str=plotForm2Par(p,fr,pd,0);
          else if(strcmp(p->fitType,"2parpoly3")==0)
            str=plotForm2ParPoly3(p,fr,pd,0);
          gnuplot_plot_equation(handle, str, "Fit");
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
