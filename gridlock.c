//definitions
#include "gridlock.h"
//functions
#include "import_data.c"
#include "print_data_info.c"
#include "generate_sums.c"
#include "print_results.c"
#include "plot_data.c"
//fitting routines
#include "1parfit.c"
#include "2parfit.c"
#include "3parfit.c"

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
  if(p->numVar==1) //see 1parfit.c
    {
      fit1Par(d,fr);
      if(strcmp(p->dataType,"chisq")==0)
        fit1ParChisqConf(fr);//generate confidence interval bounds for chisq data
    }
  else if(p->numVar==2) //see 2parfit.c
    {
      fit2Par(d,fr);
      if(strcmp(p->dataType,"chisq")==0)
        fit2ParChisqConf(fr);//generate confidence interval bounds for chisq data
    }
  else if(p->numVar==3) //see 3parfit.c
    {
      fit3Par(d,fr);
      if(strcmp(p->dataType,"chisq")==0)
        fit3ParChisqConf(fr);//generate confidence interval bounds for chisq data
    }
  
  printResults(d,p,fr); //see print_results.c
  
  if((p->plotData==1)&&(p->verbose<1))
    {
      getPlotDataNearMin(d,p,fr,pd);
      plotData(d,p,fr,pd);
    }
  
  //free structures
  free(d);
  free(p);
  free(fr);
    
  return 0; //great success
}
