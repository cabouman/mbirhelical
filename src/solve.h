#ifndef _SOLVE_H_
#define _SOLVE_H_

float  solve(
	float (*f)(float, void *),	/* pointer to function to be solved */
	float a,	/* minimum value of solution */
	float b,	/* maximum value of solution */
	float err,	/* accuarcy of solution */
	int *code,	/* error code */
	void *params	/* parameters to be handed to function */
);

float  solve_lower(
	float (*f)(),	/* pointer to function to be solved */
	float a,	/* minimum value of solution */
	float b,	/* maximum value of solution */
	float err,	/* accuarcy of solution */
	int *code,	/* error code */
	void *params	/* parameters to be handed to function */
);

#endif
