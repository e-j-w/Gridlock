//imports data from file
void importData(data * d, parameters * p)
{

  FILE *inp;
  int i;
  char str[256],str2[256],str3[256];
  
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
        	if(sscanf(str,"%s %s %Lf",str2,str3,&p->fitOpt)==3)
            {
              if(strcmp(str2,"FIT")==0)
                strcpy(p->fitType,str3);
              else
              	p->fitOpt=0.;
            }
          else if(sscanf(str,"%s %s",str2,str3)==2)
            {
              if(strcmp(str2,"FIT")==0)
                strcpy(p->fitType,str3);
              if(strcmp(str2,"UNIFORM_WEIGHT")==0)
              	{
              		p->uniWeight=1;
              		p->uniWeightVal=(long double)atof(str3);
              		p->fitOpt=0.;
              	}
            }
					else if(strcmp(str,"NONVERBOSE\n")==0)
						p->verbose=1;//only print the fit vertex data, unless an error occurs
					else if(strcmp(str,"WEIGHTED\n")==0)
						p->readWeights=1;//data has weights, in the last column
					else if(strcmp(str,"UNWEIGHTED\n")==0)
						p->readWeights=0;//data is unweighted
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
  else if(strcmp(p->fitType,"")==0)
    {
      printf("ERROR: a fit type must be specified.\nMake sure to include a line in the file with the format\n\nFIT  type\n\nwhere 'type' is a valid fit type (eg. 'par1').\n");
      exit(-1);
    }
  else
    {
      printf("ERROR: invalid fit type '%s' specified.\n",p->fitType);
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
    }
  fclose(inp);
  
  //import data from file
  if((inp=fopen(p->filename,"r"))==NULL)
    exit(-1);
  while(!(feof(inp)))//go until the end of file is reached
    {
      if(fgets(str,256,inp)!=NULL)
        {
          if(
          ((p->numVar>0)
          &&(p->readWeights==0)
          &&(sscanf(str,"%Lf %Lf %Lf %Lf %Lf %Lf",&d->x[0][d->lines],&d->x[1][d->lines],
              &d->x[2][d->lines],&d->x[3][d->lines],&d->x[4][d->lines],
              &d->x[5][d->lines])==p->numVar+1))
          ||
          ((p->numVar>0)
          &&(p->readWeights==1)
          &&(sscanf(str,"%Lf %Lf %Lf %Lf %Lf %Lf",&d->x[0][d->lines],&d->x[1][d->lines],
              &d->x[2][d->lines],&d->x[3][d->lines],&d->x[4][d->lines],
              &d->x[5][d->lines])==p->numVar+2))
          )
            {
              lineValid=1;
              for(i=0;i<p->numVar;i++)
                if(i<POWSIZE)
                  if((d->x[i][d->lines]>p->ulimit[i])||(d->x[i][d->lines]<p->llimit[i]))//check against limits
                    lineValid=0;
                    
              //deal with weights
              if(p->uniWeight==1)
              	d->x[p->numVar+1][d->lines]=p->uniWeightVal;
              else if(p->readWeights==0)
                d->x[p->numVar+1][d->lines]=1.;//set weights to 1
              if(d->x[p->numVar+1][d->lines]<=0)
                lineValid=0;//invalidate data points with bad weights (can't divide by 0 weight)
              if(lineValid==1)
              	{
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
              		
              		//got to the next data point
                	d->lines++;
                }
              else
                invalidLines++;
            }
          else if((p->numVar>0)&&(sscanf(str,"%s %Lf %Lf %Lf %Lf %Lf",str2,&d->x[0][d->lines],
              &d->x[1][d->lines],&d->x[2][d->lines],&d->x[3][d->lines],
              &d->x[4][d->lines])==p->numVar+1))
            {
              if(strcmp(str2,"UPPER_LIMITS")==0)
                {
                  for(i=0;i<p->numVar;i++)
                    if(i<POWSIZE)
                      p->ulimit[i]=d->x[i][d->lines];
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
                  for(i=0;i<p->numVar;i++)
                    if(i<POWSIZE)
                      p->llimit[i]=d->x[i][d->lines];
                  if(p->verbose<1)
                    {
                      printf("Set fit region lower limits to [");
                      for(i=0;i<p->numVar;i++)
                        printf(" %0.3LE ",p->llimit[i]);
                      printf("]\n");
                    }
                }
            }
          else if(sscanf(str,"%s %s",str2,str3)==2)
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
                  p->plotData=1;
                  strcpy(p->dataType,str3);
                  if(p->verbose<1)
                    if(strcmp(p->dataType,"chisq")==0)
                      printf("Will treat data points as chi-squared values.\n");
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
              else if(p->verbose<1)
                printf("WARNING: Improperly formatted data on line %i of the input file.\n",linenum+1);
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
      printf("ERROR: no data could be read from the input file.\n");
      if(invalidLines>0)
        printf("%i lines were skipped due to the fit region limits specified in the file.  Consider changing these limits.\n",invalidLines);
      exit(-1);
    }
  else if(p->verbose<1)
    {
      printf("Successfully read data file: %s\n%i lines of data used.\n",p->filename,d->lines);
      if(invalidLines>0)
        printf("%i lines of data skipped (outside of fit region limits).\n",invalidLines);
    }
  
}
