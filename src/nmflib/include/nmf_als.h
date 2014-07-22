#ifndef UNIVIE_NMFLIB_NMF_ALS_H
#define UNIVIE_NMFLIB_NMF_ALS_H

 /**
 * NMF_ALS -	calculates the nmf using an alternative least squares method
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
 *		The alternative least squares method 
 *		is based on following algorithm (in Matlab notation):
 *		
 *		W = rand(m,k);
 *		for iter=1:maxiter
 *			H = max(W0\A, 0);
 *			W = (H'\A')';
 *			W = max(W, 0);
 *		end
 *
 *		The matrix right division is implemented by following Lapack routines, which solve the linear system
 *
 *		W0 * x = A (H' * x = A')
 *
 *		by reducing the linear system to a triangular linear system R*x = P*Q'*A (R*x = P*Q'*A') and solving
 *		that system. P is the permutation matrix of the QR factorisation.
 *		
 *		dgeqp3	computes a QR factorisation with col. pivoting of matrix W0 (H') and stores it in W0 (H')
 *		dormqr	multiplies the matrix Q from dgeqp3 with another matrix
 *		dtrtrs	solves a triangular linear system
 *			
 * Arguments:
 *
 * a		in, 	pointer to matrix to factorise
 *
 * w0		in, 	pointer to initialised factor-matrix w0
 *		out, 	pointer to final factor-matrix w
 * h0		in, 	pointer to initialised factor-matrix h0
 *		out, 	pointer to final factor-matrix h
 * m		in, 	first dimension of a and w/w0
 *
 * n		in, 	second dimension of a and h/h0
 *
 * k		in, 	second dimension of w/w0 and first dimension of h/h0
 *
 * maxiter	in, 	maximal number of iterations to run
 *		out, 	number of iterations actually run
 * TolX		in, 	used in check for convergence, tolerance for maxchange
 *
 * TolFun	in, 	used in check for convergence, tolerance for dnorm
 *
 */


double nmf_als(double ** a, double ** w0, double ** h0, int m, int n,int k, int * maxiter, const double TolX, const double TolFun);


#endif
