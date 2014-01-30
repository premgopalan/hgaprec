#ifndef UNIVIE_NMFLIB_CHECKARGUMENTS_H
#define UNIVIE_NMFLIB_CHECKARGUMENTS_H


/**
 * checkArguments - performs basic checks on the arguments passed to nmfDriver
 *
 * Purpose:
 *		Performs sanity checks on the arguments passed to nmfDriver to check for inconsistencies and errors
 *		of the passed arguments
 *
 * Arguments:
 *
 * a		in, 	filename to load matrix a from
 *
 * k		in, 	approximation factor
 *
 * iter		in, 	maximal number of iterations to perform
 *
 * rep		in, 	number of repeated factorisations with differently initialised matrices w0 and h0
 *
 * alg		in, 	algorithm to use for the factorisation
 *
 * init		in, 	defines how to initialise the matrices w0 and h0
 *
 * min_init	in, 	defines the minimum value for initialisation
 *
 * max_init	in, 	defines the maximum value for initialisation
 *
 * w0		in, 	filename to load matrix w0 from
 *
 * h0		in, 	filename to load matrix h0 from
 *
 * w_out	in, 	filename to store final matrix w to
 *
 * h_out	in, 	filename to store final matrix h to
 *
 * TolX		in, 	tolerance value for convergence check of maxChange
 *
 * TolFun	in, 	tolerance value for convergence check of dnorm
 */

int checkArguments(const char* a, const int k, int iter, const char* w0, const char* h0, options_t * opts);


#endif
