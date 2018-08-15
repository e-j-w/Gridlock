//definitions
#include "gridlock.h"
//common functions
#include "import_data.c"
#include "print_data_info.c"
#include "generate_sums.c"
#include "plot_data.c"
//data filters
#include "lin_filter.c"
//fitting routines
#include "linfit.c"
#include "linfit_deming.c"
#include "1parfit.c"
#include "2parfit.c"
#include "3parfit.c"
#include "poly3fit.c"
#include "poly4fit.c"
#include "2parpoly3fit.c"

int main(int argc, char *argv[])
{

	//set up handler to take action upon SIGINT (CTRL-C command)
	struct sigaction sigIntHandler;
	sigIntHandler.sa_handler = sigint_cleanup;
	sigaction(SIGINT, &sigIntHandler, NULL);

	if(argc!=2)
		{
			printf("\ngridlock filename\n-----------------\n\n");
			printf("Fits the data in the specified file.\nThe file should be in plaintext, with a line specifying the fit type using the format:\n\nFIT  type\n");
			printf("\nPossible values of 'type' are:\n\nlin (line)\nlin_deming (line with errors in x)\npoly1 (1st order polynomial)\npoly2 (2nd order polynomial)\npoly3 (3rd order polynomial)\npar1 (2nd order polynomial)\npar2 (2nd order bivariate polynomial)\npar3 (2nd order trivariate polynomial)\n2parpoly3 (3rd order bivariate polynomial)\n");
			printf("\nSee the README for more details.\n");
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
  
	if(p->filter==1)
  		linearFilter(d,p);
  
	if(p->verbose<1)
		printDataInfo(d,p); //see print_data_info.c

	generateSums(d,p); //construct sums for fitting (see generate_sums.c) 
		
	//Call specific fitting routines depending on the fit type specified
	if(strcmp(p->fitType,"par1")==0) //see 1parfit.c
		fit1Par(p,d,fr,pd,1);
	else if(strcmp(p->fitType,"par2")==0) //see 2parfit.c
		fit2Par(p,d,fr,pd,1);
	else if(strcmp(p->fitType,"par3")==0) //see 3parfit.c
		fit3Par(p,d,fr,pd,1);
	else if(strcmp(p->fitType,"lin")==0)
		fitLin(p,d,fr,pd,1);
	else if(strcmp(p->fitType,"lin_deming")==0)
		{
			if(p->fitOpt==0.)//default value
				p->fitOpt=1.;
			fitLinDeming(p,d,fr,pd,1);
		}
	else if(strcmp(p->fitType,"poly3")==0)
		fitPoly3(p,d,fr,pd,1);
	else if(strcmp(p->fitType,"poly4")==0)
		fitPoly4(p,d,fr,pd,1);
	else if(strcmp(p->fitType,"2parpoly3")==0)
		fit2ParPoly3(p,d,fr,pd,1);
	
	//free structures
	free(d);
	free(p);
	free(fr);
	free(pd);
		
	return 0; //great success
}
