#ifndef UNIVIE_NMFLIB_SETDEFAULTOPTS_H
#define UNIVIE_NMFLIB_SETDEFAULTOPTS_H

/**
 * set_default_opts - Sets the default options in an options_t structure
 *
 * Purpose:
 *		If a null-pointer is passed to nmfDriver for the options-structure default values are used
 *		wich this routine implements
 * 
 * Description:
 *		The options structure contains following elements:
 *		rep		number of repetitions to run with differently initialised matrices	1
 *		init		method to use for initialising matrices					ran
 *		min_init	minimal value for random numbers in initialised matrices		0
 *		max_init	maximal value for random numbers in initialised matrices		1
 *		w_out		filename to store final matrix w in					"final_w.matrix"
 *		h_out		filename to store final matrix h in					"final_h.matrix"
 *		TolX		tolerance limit for maxChange						1E-04
 *		TolFun		tolerance for root mean square residual					1E-04
 *		nndsvd_maxiter	maximum iterations for SVD in ddsvd initialisation 			-1 -> default value set in generateMatrix
 *		nnsvd_blocksize	blocksize for SVD in ddsvd initialisation				64
 *		nndsvd_tol	tolerance for SVD in ddsvd initialisation				2E-16
 *		nndsvd_ncv	largest number of basis vectors in the Arnoldi Process			-1 -> ncv = 2 * nev (will be set in generateMatrix)
 *
 * Arguments:
 *
 * opts		in/out, 	option structure to store the default options in
 */

void set_default_opts(options_t * opts);



#endif

