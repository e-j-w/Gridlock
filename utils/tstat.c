#include "tstat.h"

//finds the approximate Student's t-score for the given input parameters
//using a lookup table
long double t_stat(int df,long double a)
{
	long double ts0[33] = {3.078,
1.886,
1.638,
1.533,
1.476,
1.44,
1.415,
1.397,
1.383,
1.372,
1.363,
1.356,
1.35,
1.345,
1.341,
1.337,
1.333,
1.33,
1.328,
1.325,
1.323,
1.321,
1.319,
1.318,
1.316,
1.315,
1.314,
1.313,
1.311,
1.31,
1.296,
1.289,
1.282
  };
  
  long double ts1[33] = {6.314,
2.92,
2.353,
2.132,
2.015,
1.943,
1.895,
1.86,
1.833,
1.812,
1.796,
1.782,
1.771,
1.761,
1.753,
1.746,
1.74,
1.734,
1.729,
1.725,
1.721,
1.717,
1.714,
1.711,
1.708,
1.706,
1.703,
1.701,
1.699,
1.697,
1.671,
1.658,
1.645
  };
  
  long double ts2[33] = {12.706,
4.303,
3.182,
2.776,
2.571,
2.447,
2.365,
2.306,
2.262,
2.228,
2.201,
2.179,
2.16,
2.145,
2.131,
2.12,
2.11,
2.101,
2.093,
2.086,
2.08,
2.074,
2.069,
2.064,
2.06,
2.056,
2.052,
2.048,
2.045,
2.042,
2,
1.98,
1.96
  };
  
  long double a_ind,df_ind;
  a_ind=-1;
  if(a==0.1)
  	a_ind=0.;
  else if(a==0.05)
  	a_ind=1.;
  else if(a==0.025)
  	a_ind=2.;
  else
  	printf("WARNING: improper value of a (%LF) passed to the t-statistic calculator.\n",a);

  if(df>0 && df<=30)
  	{
  		if(a_ind==0.)
  			return ts0[df-1];
  		else if(a_ind==1.)
  			return ts1[df-1];
  		else if(a_ind==2.)
  			return ts2[df-1];
  	}
  else if(df>30 && df<60)//interpolate
  	{
  		if(a_ind==0.)
  			return ts0[29] - (ts0[29]-ts0[30])*(df-30)/30;
  		else if(a_ind==1.)
  			return ts1[29] - (ts1[29]-ts1[30])*(df-30)/30;
  		else if(a_ind==2.)
  			return ts2[29] - (ts2[29]-ts2[30])*(df-30)/30;
  	}
  else if(df==60)
  	{
  		if(a_ind==0.)
  			return ts0[30];
  		else if(a_ind==1.)
  			return ts1[30];
  		else if(a_ind==2.)
  			return ts2[30];
  	}
  else if(df>60 && df<120)//interpolate
  	{
  		if(a_ind==0.)
  			return ts0[30] - (ts0[30]-ts0[31])*(df-60)/60;
  		else if(a_ind==1.)
  			return ts1[30] - (ts1[30]-ts1[31])*(df-60)/60;
  		else if(a_ind==2.)
  			return ts2[30] - (ts2[30]-ts2[31])*(df-60)/60;
  	}
  else if(df==120)
  	{
  		if(a_ind==0.)
  			return ts0[31];
  		else if(a_ind==1.)
  			return ts1[31];
  		else if(a_ind==2.)
  			return ts2[31];
  	}
  else if(df>120)//use value at infinity
  	{
  		if(a_ind==0.)
  			return ts0[32];
  		else if(a_ind==1.)
  			return ts1[32];
  		else if(a_ind==2.)
  			return ts2[32];
  	}
  
  printf("WARNING: improper parameters (df=%i,a=%LF) passed to the t-statistic calculator.  Returning 0...\n",df,a);
	return 0.;
}
