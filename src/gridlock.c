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
#include "poly2fit.c"
#include "poly2root0fit.c"
#include "2parpoly2fit.c"
#include "3parpoly2fit.c"
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
			printf("usage: gridlock filename\n\n");
			printf("Fits the data in the plaintext file specified by 'filename'.\nThe fit type and data should be specified in the file using the format:\n\nFIT  type\nVariableValue1  DataValue1\nVariableValue2  DataValue2\n...             ...\n");
			printf("\nPossible values of 'type' are:\nlin (linear / 1st order polynomial)\nlin_deming (linear with errors in x)\npoly2 (2nd order polynomial)\npoly3 (3rd order polynomial)\npoly4 (4th order polynomial)\npoly2root0 (2nd order polynomial, 0th order term fixed to zero)\n2parpoly2 (2nd order bivariate polynomial)\n2parpoly3 (3rd order bivariate polynomial)\n3parpoly2 (2nd order trivariate polynomial)\n");
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
			printf("ERROR: the number of free parameters (NUM_PAR) must be 3 or less, and cannot be negative.\n");
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
	if(strcmp(p->fitType,"poly2")==0) //see poly2fit.c
		fit1Par(p,d,fr,pd,1);
	else if(strcmp(p->fitType,"poly2root0")==0) //see poly2root0fit.c
		fitPoly2Root0(p,d,fr,pd,1);
	else if(strcmp(p->fitType,"2parpoly2")==0) //see 2parpoly2fit.c
		fit2Par(p,d,fr,pd,1);
	else if(strcmp(p->fitType,"3parpoly2")==0) //see 3parpoly2fit.c
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
