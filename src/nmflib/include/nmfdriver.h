#ifndef UNIVIE_NMFLIB_NMFDRIVER_H
#define UNIVIE_NMFLIB_NMFDRIVER_H

#include "common.h"

/**
 * nmfDriver - performs a whole factorisation process including necessary steps before and after the factorisation
 * 
 * Purpose:
 *		This routine calculates a non negative matrix factorisation of a m by n matrix A
 *
 *		A = W * H
 *
 *		where A, W and H are non-negative matrices.
 *		W is a m by k matrix
 *		H is a k by n matrix
 *		k is the approximation factor
 *		
 * Description:
 *		This routine loads a matrix A (and in case they are specified matrices w0 and h0) from a file and
 *		performs multiple iterations of one of the implemented nmf-algorithms, until an iteration or tolerance
 *		limit is reached.
 *		The initial factor matrices W0 and H0 can be loaded from a file or initialised randomly.
 *		It is possible to run several complete runs, each with different inital matrices w0 and h0. For every run
 *		after the first the matrices then get randomly initialised.
 *		
 *		The factorisation results are normalised and rows of H and columns of W are resorted. The order is
 *		determined by the row sums of elementwise squared matrix H.
 *		
 *		Factorisation results of different runs are compared by the root mean square residual and the best result
 *		of all runs is stored in a file.
 *		Caution: The root mean square residual in this routine is the frobenius norm divided by sqrt(m*n)
 *
 *		Implemented algorithms to actually compute the factorisation are (for details see each algorithm):
 *
 *		nmf_mu		multiplicative update approach
 *		nmf_als		alternating least squares approach
 *		nmf_neals	normal equation alternating least squares approach
 *		nmf_alspg	alternating least squares approach using a projected gradient method
 *		nmf_pg		direct projected gradient approach
 *
 *
 * Usage:
 *		For calling nmfDriver at least arguments a, k and iter have to be passed.
 *		Passing a NULL-pointer to w0 or h0 will let w0 and h0 be initialised according to opts->init
 *		Passing a NULL-pointer to opts will set default options
 *
 *		The datatype "options_t" has following structure:
 *
 * rep		in, 	number of repeated factorisations with differently initialised matrices w0 and h0
 *
 *
 * init		in, 	defines how to initialise the matrices w0 and h0
 *
 * min_init	in, 	defines the minimum value for initialisation
 *
 * max_init	in, 	defines the maximum value for initialisation
 *
 * w_out	in, 	filename to store final matrix w to
 *
 * h_out	in, 	filename to store final matrix h to
 *
 * TolX		in, 	tolerance value for convergence check of maxChange
 *			Every iteration checks for the change in each factor matrix
 * TolFun	in, 	tolerance value for convergence check of root mean square residual
 *
 *
 * Arguments:
 *
 * a		in, 	filename to load matrix a from
 *
 * k		in, 	approximation factor which limits the matrices dimensions
 *
 * iter		in, 	maximal number of iterations to perform
 *
 * w0		in, 	filename to load matrix w0 from, when NULL matrix is initialised using method in "opts->init"
 *
 * h0		in, 	filename to load matrix h0 from, when NULL matrix is initialised using method in "opts->init"
 *
 * alg		in,	algorithm to use for the factorisation steps
 *
 * opts		in,	options_t structure, which defines all additional options
 *			whenn NULL, defaultvalues are used
 */

extern "C" {
int nmfDriver(const char* a, int k, int iter, const char* w0, const char* h0, alg_t alg, options_t * opts);
}




#endif
