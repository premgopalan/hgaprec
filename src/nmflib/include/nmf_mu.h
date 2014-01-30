#ifndef UNIVIE_NMFLIB_NMF_MU_H
#define UNIVIE_NMFLIB_NMF_MU_H

/**
 * NMF_MU -	calculates the nmf using a multiplicative update method
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
 *		The multiplicative update method as described in Berry (2006) (see references [1])
 *		is based on following algorithm (in Matlab notation):
 *
 *		W = rand(m, k);
 *		H = rand(k, n);
 *		for iter = 1:maxiter
 *			H = H .* (W' * A) ./ (W' * W * H + 10E-09)
 *			W = W .* (A * H') ./ (W * H * H' + 10E-09)
 *		end
 *
 *		This routine uses the BLAS "dgemm" routine for matrix-matrix-multiplications
 *
 * Arguments:
 *
 * a		in, 	pointer to matrix to factorise
 *
 * w0		in, 	pointer to initialised factor-matrix w0
 *		out, 	pointer to final factor-matrix w
 *
 * h0		in, 	pointer to initialised factor-matrix h0
 *		out, 	pointer to final factor-matrix h
 *
 * m		in, 	first dimension of a and w/w0
 *
 * n		in, 	second dimension of a and h/h0
 *
 * k		in, 	second dimension of w/w0 and first dimension of h/h0
 *
 * maxiter	in, 	maximal number of iterations to run
 *		out, 	number of iterations actually run
 *
 * TolX		in, 	used in check for convergence, tolerance for maxchange of matrices w and h
 *
 * TolFun	in, 	used in check for convergence, tolerance for root mean square residual
 */


double nmf_mu(double * a, double ** w0, double ** h0, int m, int n, \
		      int k, int * maxiter, const double TolX, const double TolFun);



#endif
