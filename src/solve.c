#include "solve.h"

float  solve(
	float (*f)(float, void *), /* pointer to function to be solved */
	float a,      /* minimum value of solution */
	float b,      /* maximum value of solution */
	float err,    /* accuarcy of solution */
	int *code,     /* error code */
	void *params   /* parameters to be handed to function */
)
/* Solves equation (*f)(x,void *params) = 0 on x in [a,b]. Uses half interval method.*/
/* Requires that (*f)(a) and (*f)(b) have opposite signs.		*/
/* Returns code=0 if signs are opposite.				*/
/* Returns code=1 if signs are both positive. 				*/
/* Returns code=1 if signs are both negative. 				*/
{
	int     signa,signb,signc;
	float  fa,fb,fc,c,signaling_nan();
	float  dist;

	fa = (*f)(a,params);  signa = fa>0;
	fb = (*f)(b,params);  signb = fb>0;

	/* check starting conditions */
	if( signa==signb ) {
		if(signa==1) *code = 1;
		else *code = -1;
		return(0.0);
	}
	else *code = 0;

	/* half interval search */
	if( (dist=b-a)<0 ) dist = -dist;
	while(dist>err) {
		c = (b+a)/2;
		fc = (*f)(c,params);  signc = fc>0;
		if(signa == signc) { a = c; fa = fc; }
		else { b = c; fb = fc; }
		if( (dist=b-a)<0 ) dist = -dist;
	}

	/* linear interpolation */
	if( (fb-fa)==0 ) return(a);
	else {
		c = (a*fb - b*fa)/(fb-fa);
		return(c);
	}
}


/* solve_lower returns lowest estimate within tolerance */
/* while solve linearly interpolates */

float  solve_lower(
	float (*f)(), /* pointer to function to be solved */
	float a,      /* minimum value of solution */
	float b,      /* maximum value of solution */
	float err,    /* accuarcy of solution */
	int *code,     /* error code */
	void *params   /* parameters to be handed to function */
)
/* Solves equation (*f)(x,void *params) = 0 on x in [a,b]. Uses half interval method.*/
/* Requires that (*f)(a) and (*f)(b) have opposite signs.		*/
/* Returns code=0 if signs are opposite.				*/
/* Returns code=1 if signs are both positive. 				*/
/* Returns code=1 if signs are both negative. 				*/
{
	int     signa,signb,signc;
	float  fa,fb,fc,c,signaling_nan();
	float  dist;

	fa = (*f)(a,params);  signa = fa>0;
	fb = (*f)(b,params);  signb = fb>0;

	/* check starting conditions */
	if( signa==signb ) {
		if(signa==1) *code = 1;
		else *code = -1;
		return(0.0);
	}
	else *code = 0;

	/* half interval search */
	if( (dist=b-a)<0 ) dist = -dist;
	while(dist>err) {
		c = (b+a)/2;
		fc = (*f)(c,params);  signc = fc>0;
		if(signa == signc) { a = c; fa = fc; }
		else { b = c; fb = fc; }
		if( (dist=b-a)<0 ) dist = -dist;
	}

	/* linear interpolation */
	/*if( (fb-fa)==0 ) return(a);
	  else {
	  c = (a*fb - b*fa)/(fb-fa);
	  return(c);
	  } 
	  */
	return(a);
}
