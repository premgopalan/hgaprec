#ifndef UNIVIE_NMFLIB_NMF_PG_H
#define UNIVIE_NMFLIB_NMF_PG_H


/**
 * NMF_PG -	calculates the nmf using a direct projected gradient approach
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
 *		This routine implements an direct projected gradients gradients approach.
 *		Its algorithm is explained in Lin (2007) (see references [2])
 *
 *		The general gradient descent approach is shown in following code:
 *
 *		W = rand(m,k);
 *		H = rand(k,n);
 *		for iter=1:maxiter
 *			H = H - epsH * gradH;
 *			W = W - epsW * gradW;
 *		end
 *		
 *		The gradients are calculated as follows:
 *
 *		gradW = W * (H*H') - A*H';
 *		gradH = (W'*W)*H - W'*V;
 *
 *		This routine uses following BLAS/LAPACK routines:
 *		
 *		dgemm	matrix-matrix multiplications to calculate gradients
 *		dlange	used to compute the frobenius norm
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
 * tol		in, 	used in check for convergence, relative stopping criterion
 */

double nmf_pg(double ** a, double ** w0, double **h0, int m, int n, int k, int * maxiter, const double tol);


#endif
