//imports data from file
void importData(data * d, parameters * p)
{

  FILE *inp;
  int i,j;
  char str[256],str2[256],str3[256];
  long double val;
  
  //initialize values
  int invalidLines=0;
  int lineValid;
  int linenum=0;
  p->plotData=0;
  for(i=0;i<POWSIZE;i++)
    {
      p->llimit[i]=-1*BIG_NUMBER;
      p->ulimit[i]=BIG_NUMBER;
      d->max_x[i]=-1*BIG_NUMBER;
      d->min_x[i]=BIG_NUMBER;
    }
  p->dllimit=-1*BIG_NUMBER;
  p->dulimit=BIG_NUMBER;
  p->numCIEvalPts=0;
  d->max_m=-1*BIG_NUMBER;
  d->min_m=BIG_NUMBER;
    
  if((inp=fopen(p->filename,"r"))==NULL)
    {
      printf("\nERROR: input file %s can not be opened.\n",p->filename);
      exit(-1);
    }

  //read the number of parameters from the file and set verbosity of output
  while(!(feof(inp)))//go until the end of file is reached
    {
      if(fgets(str,256,inp)!=NULL)
        {
        	if(sscanf(str,"%s %s %Lf",str2,str3,&val)==3)
            {
              if(strcmp(str2,"FIT")==0){
                strcpy(p->fitType,str3);
                p->fitOpt = val;
              }else if(strcmp(str2,"SLICE_PAR")==0){
                if(strcmp(str3,"x")==0){
                  p->ignorePar[0]=2;
                  p->sliceVal[0] = val;
                }else if(strcmp(str3,"y")==0){
                  p->ignorePar[1]=2;
                  p->sliceVal[1] = val;
                }else if(strcmp(str3,"z")==0){
                  p->ignorePar[2]=2;
                  p->sliceVal[2] = val;
                }
              }
            }
          else if(sscanf(str,"%s %s",str2,str3)==2)
            {
              if(strcmp(str2,"FIT")==0)
                strcpy(p->fitType,str3);
              else if(strcmp(str2,"UNIFORM_WEIGHT")==0)
              	{
              		p->uniWeight=1;
              		p->uniWeightVal=(long double)atof(str3);
              		p->fitOpt=0.;
              	}
              else if(strcmp(str2,"LINEAR_FILTER")==0)
              	{
              		p->filter=1;//use linear filter on data
              		p->filterSigma=atof(str3);
              	}
              else if(strcmp(str2,"EVAL_CI")==0)
              	{
              		p->CIEvalPts[p->numCIEvalPts]=(long double)atof(str3);
                  p->numCIEvalPts++;
              	}
              else if(strcmp(str2,"IGNORE_PAR")==0)
              	{
                  if(strcmp(str3,"x")==0)
              		  p->ignorePar[0]=1;
                  else if(strcmp(str3,"y")==0)
                    p->ignorePar[1]=1;
                  else if(strcmp(str3,"z")==0)
                    p->ignorePar[2]=1;
              	}
            }
					else if(strcmp(str,"PARAMETERS\n")==0)
						p->verbose=1;//only print the fit vertex data, unless an error occurs
					else if(strcmp(str,"COEFFICIENTS\n")==0)
						p->verbose=2;//only print the fit coefficients, unless an error occurs
					else if((strcmp(str,"WEIGHTED\n")==0)||(strcmp(str,"WEIGHT\n")==0)||(strcmp(str,"WEIGHTS\n")==0))
						p->readWeights=1;//data has weights, in the last column
					else if(strcmp(str,"UNWEIGHTED\n")==0)
						p->readWeights=0;//data is unweighted
          else if(strcmp(str,"ZEROX\n")==0)
						p->forceZeroX=1;//force x to zero
          else if(strcmp(str,"ZEROY\n")==0)
						p->forceZeroY=1;//force y to zero
          else if(strcmp(str,"FIND_MIN_GRID_POINT_FROM_FIT\n")==0)
						p->findMinGridPoint=1;//find the grid point corresponding to the smallest value of the fit function
          else if(strcmp(str,"FIND_MAX_GRID_POINT_FROM_FIT\n")==0)
						p->findMaxGridPoint=1;//find the grid point corresponding to the smallest value of the fit function
        }
    }
  //check the fit type
  if(strcmp(p->fitType,"poly1")==0)
  	strcpy(p->fitType,"lin");
  else if(strcmp(p->fitType,"poly2")==0)
  	strcpy(p->fitType,"par1");
  if(strcmp(p->fitType,"par1")==0)
    p->numVar=1;
  else if(strcmp(p->fitType,"par2")==0)
    p->numVar=2;
  else if(strcmp(p->fitType,"par3")==0)
    p->numVar=3;
  else if(strcmp(p->fitType,"lin")==0)
  	{
    	p->numVar=1;
    	p->plotCI=1;
    }
  else if(strcmp(p->fitType,"lin_deming")==0)
    p->numVar=1;
  else if(strcmp(p->fitType,"poly3")==0)
    p->numVar=1;
  else if(strcmp(p->fitType,"poly4")==0)
    p->numVar=1;
  else if(strcmp(p->fitType,"2parpoly3")==0)
    p->numVar=2;
  else if(strcmp(p->fitType,"")==0)
    {
      printf("ERROR: a fit type must be specified.\nMake sure to include a line in the file with the format\n\nFIT  type\n\nwhere 'type' is a valid fit type (eg. 'par1').\n");
      printf("\nValid fit types are:\n\nlin (line)\nlin_deming (line with errors in x)\npoly1 (1st order polynomial)\n");
      printf("poly2 (2nd order polynomial)\npoly3 (3rd order polynomial)\npoly4 (4th order polynomial)\npar1 (2nd order polynomial)\n");
      printf("par2 (2nd order bivariate polynomial)\npar3 (2nd order trivariate polynomial)\n2parpoly3 (3rd order bivariate polynomial)\n");
      exit(-1);
    }
  else
    {
      printf("ERROR: invalid fit type '%s' specified.\n",p->fitType);
      printf("\nValid fit types are:\n\nlin (line)\nlin_deming (line with errors in x)\npoly1 (1st order polynomial)\n");
      printf("poly2 (2nd order polynomial)\npoly3 (3rd order polynomial)\npoly4 (4th order polynomial)\npar1 (2nd order polynomial)\n");
      printf("par2 (2nd order bivariate polynomial)\npar3 (2nd order trivariate polynomial)\n2parpoly3 (3rd order bivariate polynomial)\n");
      exit(-1);
    }
    
  if(p->verbose<1)
    {
      if(strcmp(p->fitType,"par1")==0)
        printf("Will fit a parabola.\n");
      if(strcmp(p->fitType,"par2")==0)
        printf("Will fit a paraboloid with %i free parameters.\n",p->numVar);
      if(strcmp(p->fitType,"par3")==0)
        printf("Will fit a paraboloid with %i free parameters.\n",p->numVar);
      if(p->uniWeight==1)
      	printf("Uniform weights of value %0.3Lf will be taken.\n",p->uniWeightVal);
      else if(p->readWeights==0)
        printf("No weights will be taken for data points.\n");
      else if(p->readWeights==1)
        printf("Weights for data points will be taken from the last column of the data file.\n");
      if(p->ignorePar[0]==1)
        printf("Will ignore data corresponding to the x variable.\n");
      if(p->ignorePar[1]==1)
        printf("Will ignore data corresponding to the y variable.\n");
      if(p->ignorePar[2]==1)
        printf("Will ignore data corresponding to the z variable.\n");
      if(p->ignorePar[0]==2)
        printf("Will slice data at x = %Lf.\n",p->sliceVal[0]);
      if(p->ignorePar[1]==2)
        printf("Will slice data at y = %Lf.\n",p->sliceVal[1]);
      if(p->ignorePar[2]==2)
        printf("Will slice data at z = %Lf.\n",p->sliceVal[2]);
    }
  fclose(inp);
  
  //by default, use the appropriate 1-sigma confidence level
  strcpy(p->ciSigmaDesc,"1-sigma (68.3%)");
  if(p->numVar==1)
  	p->ciDelta=1.00;
  else if(p->numVar==2)
  	p->ciDelta=2.30;
  else if(p->numVar==3)
  	p->ciDelta=3.53;
  else
  	p->ciDelta=0.00;
  
  
  //setup data for if parameters are ignored/sliced
  int numIgnoredPar = 0;
  for(i=0;i<POWSIZE;i++){
    if(p->ignorePar[i]>=1){
      numIgnoredPar++;
    }
  }
  
  //import data from file
  int numCols;
  if((inp=fopen(p->filename,"r"))==NULL)
    exit(-1);
  while(!(feof(inp)))//go until the end of file is reached
    {
      if(fgets(str,256,inp)!=NULL)
        {
          numCols = sscanf(str,"%Lf %Lf %Lf %Lf %Lf %Lf",&d->x[0][d->lines],&d->x[1][d->lines],&d->x[2][d->lines],&d->x[3][d->lines],&d->x[4][d->lines],&d->x[5][d->lines]);
          if( ((p->numVar>0)&&(p->readWeights==0)&&(numCols==p->numVar+1+numIgnoredPar)) || ((p->numVar>0)&&(p->readWeights==1)&&(numCols==p->numVar+2+numIgnoredPar)) )
            {
              lineValid=1;

              //handle validity of sliced data
              for(i=POWSIZE-1;i>=0;i--){
                if(p->ignorePar[i]==2){
                  if(numCols > i+1){
                    if(d->x[i][d->lines] != p->sliceVal[i]){
                      lineValid = 0;
                      break;
                    }
                  }
                }
              }
              
              if(lineValid == 1){

                //handle ignored/sliced variables by reshuffling data
                for(i=POWSIZE-1;i>=0;i--){
                  if(p->ignorePar[i]>=1){
                    if(numCols > i+1){
                      for(j=i;j<p->numVar+numIgnoredPar;j++)
                        d->x[j][d->lines]=d->x[j+1][d->lines];
                      if(p->readWeights==1)
                        d->x[p->numVar+numIgnoredPar][d->lines]=d->x[p->numVar+1+numIgnoredPar][d->lines];
                    }
                  }
                }

                //check variable and data values for NaN
                for(i=0;i<p->numVar+2;i++)
                  if(i<POWSIZE)
                    if(d->x[i][d->lines]!=d->x[i][d->lines]){
                      lineValid=0;
                      break;
                    }
                      

                //check variable values against limits
                for(i=0;i<p->numVar;i++)
                  if(i<POWSIZE)
                    if((d->x[i][d->lines]>p->ulimit[i])||(d->x[i][d->lines]<p->llimit[i])){
                      lineValid=0;
                      break;
                    }
                      
                
                //check data values against limits
                if((d->x[p->numVar][d->lines]>p->dulimit)||(d->x[p->numVar][d->lines]<p->dllimit)){
                  lineValid=0;
                }

              }
              
              if(lineValid==1)
              	{

                  //deal with weights
                  if(p->uniWeight==1)
                    d->x[p->numVar+1][d->lines]=p->uniWeightVal;
                  else if(p->readWeights==0)
                    d->x[p->numVar+1][d->lines]=1.;//set weights to 1
                  if(d->x[p->numVar+1][d->lines]<=0)
                    lineValid=0;//invalidate data points with bad weights (can't divide by 0 weight)

              		//determine maximum and minimum values
              		if(d->x[p->numVar][d->lines] > d->max_m)
              			d->max_m=d->x[p->numVar][d->lines];
              		if(d->x[p->numVar][d->lines] < d->min_m)
              			d->min_m=d->x[p->numVar][d->lines];
              		for(i=0;i<p->numVar;i++)
              			{
              				if(d->x[i][d->lines] > d->max_x[i])
				          			d->max_x[i]=d->x[i][d->lines];
				          		if(d->x[i][d->lines] < d->min_x[i])
				          			d->min_x[i]=d->x[i][d->lines];
              			}
              		
              		//go to the next data point
                	d->lines++;
                }
              else
                invalidLines++;
            }
          else if(sscanf(str,"%s %s",str2,str3)>=2)
            {
              numCols = sscanf(str,"%s %Lf %Lf %Lf %Lf %Lf",str2,&d->x[0][d->lines],&d->x[1][d->lines],&d->x[2][d->lines],&d->x[3][d->lines],&d->x[4][d->lines]);
              if((p->numVar>0)&&(numCols==p->numVar+1+numIgnoredPar))
                {
                  if(strcmp(str2,"UPPER_LIMITS")==0)
                    {
                      for(i=0;i<p->numVar+numIgnoredPar;i++)
                        if(i<POWSIZE)
                          p->ulimit[i]=d->x[i][d->lines];
                      
                      //reshuffle limit if parameters ignored
                      for(i=POWSIZE-1;i>=0;i--){
                        if(p->ignorePar[i]>=1){
                          if(numCols > i+1){
                            for(j=i;j<p->numVar+numIgnoredPar;j++)
                              if(j<POWSIZE-1)
                                p->ulimit[j]=p->ulimit[j+1];
                          }
                        }
                      }

                      if(p->verbose<1)
                        {
                          printf("Set fit region upper limits to [");
                          for(i=0;i<p->numVar;i++)
                            printf(" %0.3LE ",p->ulimit[i]);
                          printf("]\n");
                        }
                    }
                  if(strcmp(str2,"LOWER_LIMITS")==0)
                    { 
                      for(i=0;i<p->numVar+numIgnoredPar;i++)
                        if(i<POWSIZE)
                          p->llimit[i]=d->x[i][d->lines];

                      //reshuffle limit if parameters ignored
                      for(i=POWSIZE-1;i>=0;i--){
                        if(p->ignorePar[i]>=1){
                          if(numCols > i+1){
                            for(j=i;j<p->numVar+numIgnoredPar;j++)
                              if(j<POWSIZE-1)
                                p->llimit[j]=p->llimit[j+1];
                          }
                        }
                      }

                      if(p->verbose<1)
                        {
                          printf("Set fit region lower limits to [");
                          for(i=0;i<p->numVar;i++)
                            printf(" %0.3LE ",p->llimit[i]);
                          printf("]\n");
                        }
                    }
                }
              if(sscanf(str,"%s %s",str2,str3)==2)
                {
                  
                  if(strcmp(str2,"PLOT")==0)
                    {
                      p->plotData=1;
                      strcpy(p->plotMode,str3);
                      if(p->verbose<1)
                        printf("Will plot data using mode: %s\n",p->plotMode);
                    }
                  if(strcmp(str2,"DATA_TYPE")==0)
                    {
                      strcpy(p->dataType,str3);
                      if(p->verbose<1)
                        if(strcmp(p->dataType,"chisq")==0)
                          printf("Will treat data points as chi-squared values.\n");
                    }
                  if(strcmp(str2,"DATA_UPPER_LIMIT")==0)
                    {
                      if(sscanf(str3,"%Lf",&p->dulimit))
                        printf("Set data upper limit to: %0.3LE\n",p->dulimit);
                      else
                        {
                          printf("ERROR: could not properly set data upper limit (DATA_UPPER_LIMIT option).\n");
                          exit(-1);
                        }
                    }
                  if(strcmp(str2,"DATA_LOWER_LIMIT")==0)
                    {
                      if(sscanf(str3,"%Lf",&p->dllimit))
                        printf("Set data lower limit to: %0.3LE\n",p->dllimit);
                      else
                        {
                          printf("ERROR: could not properly set data lower limit (DATA_LOWER_LIMIT option).\n");
                          exit(-1);
                        }
                    }
                  if(strcmp(str2,"SET_CI_DELTA")==0)
                    {
                      if(sscanf(str3,"%Lf",&p->ciDelta))
                        {
                          printf("Set confidence interval delta value to: %0.3LE\n",p->ciDelta);
                          sprintf(p->ciSigmaDesc,"custom (delta=%Lf)",p->ciDelta);//indicate custom confidence interval
                        }
                      else
                        {
                          printf("ERROR: could not properly set confidence interval delta value (SET_CI_DELTA option).\n");
                          exit(-1);
                        }
                        
                    }
                  if(strcmp(str2,"SET_CI_SIGMA")==0)
                    {
                      
                      if(strcmp(str3,"1")==0)
                        {
                          if(p->numVar==1)
                            p->ciDelta=1.00;
                          else if(p->numVar==2)
                            p->ciDelta=2.30;
                          else if(p->numVar==3)
                            p->ciDelta=3.53;
                          else
                            p->ciDelta=0.00;
                          printf("Set confidence interval to 1-sigma (68.3%%), delta value: %0.3LE\n",p->ciDelta);
                          strcpy(p->ciSigmaDesc,"1-sigma (68.3%)");
                        }
                      else if(strcmp(str3,"2")==0)
                        {
                          if(p->numVar==1)
                            p->ciDelta=4.00;
                          else if(p->numVar==2)
                            p->ciDelta=6.17;
                          else if(p->numVar==3)
                            p->ciDelta=8.02;
                          else
                            p->ciDelta=0.00;
                          printf("Set confidence interval to 2-sigma (95.4%%), delta value: %0.3LE\n",p->ciDelta);
                          strcpy(p->ciSigmaDesc,"2-sigma (95.4%)");
                        }
                      else if(strcmp(str3,"3")==0)
                        {
                          if(p->numVar==1)
                            p->ciDelta=9.00;
                          else if(p->numVar==2)
                            p->ciDelta=11.8;
                          else if(p->numVar==3)
                            p->ciDelta=14.2;
                          else
                            p->ciDelta=0.00;
                          printf("Set confidence interval to 3-sigma (99.73%%), delta value: %0.3LE\n",p->ciDelta);
                          strcpy(p->ciSigmaDesc,"3-sigma (99.73%)");
                        }
                      else if(strcmp(str3,"90%")==0)
                        {
                          if(p->numVar==1)
                            p->ciDelta=2.71;
                          else if(p->numVar==2)
                            p->ciDelta=4.61;
                          else if(p->numVar==3)
                            p->ciDelta=6.25;
                          else
                            p->ciDelta=0.00;
                          printf("Set confidence interval to 90%%, delta value: %0.3LE\n",p->ciDelta);
                          strcpy(p->ciSigmaDesc,"90%");
                        }
                      else
                        {
                          printf("ERROR: Invalid parameter for SET_CI_SIGMA: %s\nValid parameters: 1, 2, 3, 90%%\n",str3);
                          exit(-1);
                        }
                        
                    }
                }
              
            }
          else
            {
              if(strcmp(str,"PLOT\n")==0)
                {
                  p->plotData=1;
                  if(p->verbose<1)
                    printf("Will plot data.\n");
                }
              else if((strcmp(str,"PARAMETERS\n")!=0)&&(strcmp(str,"COEFFICIENTS\n")!=0)&&(strcmp(str,"WEIGHTED\n")!=0)&&
                      (strcmp(str,"WEIGHT\n")!=0)&&(strcmp(str,"WEIGHTS\n")!=0)&&(strcmp(str,"UNWEIGHTED\n")!=0)&&
                      (strcmp(str,"ZEROX\n")!=0)&&(strcmp(str,"ZEROY\n")!=0)&&(strcmp(str,"FIND_MIN_GRID_POINT_FROM_FIT\n")!=0)&&(strcmp(str,"FIND_MAX_GRID_POINT_FROM_FIT\n")!=0))
                if(p->verbose<1)
                  printf("WARNING: Improperly formatted data on line %i of the input file.\nLine content: %s",linenum+1,str);
            }
          if((p->numVar==1)&&(sscanf(str,"%s %s",str2,str3)==2))//workaround to allow plotting lines to be read when using 1 free parameter
            {
              if(strcmp(str2,"PLOT")==0)
                {
                  p->plotData=1;
                  strcpy(p->plotMode,str3);
                  if(p->verbose<1)
                    printf("Will plot data using mode: %s\n",p->plotMode);
                }
            }
          linenum++;
          if(linenum>MAXFILELENGTH)
          	{
          		printf("ERROR: Number of data points in input file exceeds MAXFILELENGTH (%i).\n",MAXFILELENGTH);
          		printf("Please reduce the number of data points used, or modify MAXFILELENGTH in gridlock.h and recompile.\n");
          		exit(-1);
          	}
        }
    }
  fclose(inp);

  if(d->lines<1)
    {
      sprintf(str,"specified fit type '%s' requires data using %i parameter(s).",p->fitType,p->numVar);
      printf("ERROR: no data could be read from the input file.  If data is present, check that it uses the correct number of parameters - %s",str);
      if(numIgnoredPar > 0){
        sprintf(str,"  Note that you are ignoring and/or slicing on %i parameter(s) in the data (IGNORE_PAR and/or SLICE_PAR commands).",numIgnoredPar);
        printf("%s",str);
      }
      printf("\n");
      if(invalidLines>0)
        printf("%i lines were skipped due to the fit region limits specified in the file.  Consider changing these limits.\n",invalidLines);
      exit(-1);
    }
  else if(p->verbose<1)
    {
      printf("Successfully read data file: %s\n%i line(s) of data used.\n",p->filename,d->lines);
      if(invalidLines>0)
        printf("%i line(s) of data skipped (outside of fit region limits).\n",invalidLines);
    }
  
}
