//definitions
#include "gridlock.h"
//fitting routines
#include "linfit.c"
#include "linfit_deming.c"
#include "1parfit.c"
#include "2parfit.c"
#include "3parfit.c"
#include "poly3fit.c"
#include "2parpoly3fit.c"
//functions
#include "import_data.c"
#include "print_data_info.c"
#include "generate_sums.c"
#include "plot_data.c"

int main(int argc, char *argv[])
{

  //set up handler to take action upon SIGINT (CTRL-C command)
  struct sigaction sigIntHandler;
  sigIntHandler.sa_handler = sigint_cleanup;
  sigaction(SIGINT, &sigIntHandler, NULL);

  if(argc!=2)
    {
      printf("\ngridlock filename\n-----------------\n\n");
      printf("Fits the data in the specified file.  The file should be in plaintext, with a line specifying the fit type using the format:\n\nFIT  type\n\nPossible values of 'type' are listed in the README.\n\n");
      exit(-1);
    }

  //allocate structures
  parameters *p=(parameters*)calloc(1,sizeof(parameters));
  data *d=(data*)calloc(1,sizeof(data));
  fit_results *fr=(fit_results*)calloc(1,sizeof(fit_results));
  plot_data *pd=(plot_data*)calloc(1,sizeof(plot_data));

  strcpy(p->filename,argv[1]);
  importData(d,p); //see import_data.c
  
  if((p->numVar<1)||(p->numVar>3))
    {
      printf("ERROR: the number of free parameters (NUM_PAR) must be 3 or less, and cannot be negative.\nAborting...\n");
      exit(-1);
    }
  if(p->numVar>(POWSIZE-2))
    {
      printf("ERROR: the number of free parameters is greater than POWSIZE - 2 (%i).\nPlease edit the value of POWSIZE in gridlock.h and recompile.\n",POWSIZE-2);
      exit(-1);
    }
  
  if(p->verbose<1)
    printDataInfo(d,p); //see print_data_info.c

  generateSums(d,p); //construct sums for fitting (see generate_sums.c) 
  
  //Call specific fitting routines depending on 
  //the number of free parameters and other settings.
  if(strcmp(p->fitType,"par1")==0) //see 1parfit.c
    {
      fit1Par(d,fr);
      if(strcmp(p->dataType,"chisq")==0)
        fit1ParChisqConf(p,fr);//generate confidence interval bounds for chisq data
      print1Par(d,p,fr);
    }
  else if(strcmp(p->fitType,"par2")==0) //see 2parfit.c
    {
      fit2Par(d,fr);
      if(strcmp(p->dataType,"chisq")==0)
        fit2ParChisqConf(p,fr);//generate confidence interval bounds for chisq data
      print2Par(d,p,fr);
    }
  else if(strcmp(p->fitType,"par3")==0) //see 3parfit.c
    {
      fit3Par(d,fr);
      if(strcmp(p->dataType,"chisq")==0)
        fit3ParChisqConf(p,fr);//generate confidence interval bounds for chisq data
      print3Par(d,p,fr);
    }
	else if(strcmp(p->fitType,"lin")==0)
		fitLin(p,d,fr,1);
	else if(strcmp(p->fitType,"lin_deming")==0)
		{
			if(p->fitOpt==0.)//default value
				p->fitOpt=1.;
			fitLinDeming(p,d,fr,p->fitOpt);
			printLinDeming(d,p,fr);
		}
	else if(strcmp(p->fitType,"poly3")==0)
		fitPoly3(p,d,fr,1);      
	else if(strcmp(p->fitType,"2parpoly3")==0)
		fit2ParPoly3(p,d,fr,1);
  
  if((p->plotData==1)&&(p->verbose<1))
    {
      preparePlotData(d,p,fr,pd);
      plotData(p,fr,pd);
    }
  
  //free structures
  free(d);
  free(p);
  free(fr);
  free(pd);
    
  return 0; //great success
}
