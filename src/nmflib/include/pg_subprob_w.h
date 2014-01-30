#ifndef UNIVIE_NMFLIB_PG_SUBPROB_W_H
#define UNIVIE_NMFLIB_PG_SUBPROB_W_H

/**
 * PG_SUBPROB_W	calculates the projected gradient step of factor matrix W used in nmf_pg and nmf_alspg
 * 
 * Purpose:
 *		This routine calculates the projected gradient step of factor matrix W used in nmf_pg and nmf_alspg
 *		It is split in two routines for W and H respectively to avoid unecessary transposition of input and output matrices
 *
 * Description:
 *		The general gradient descent approach is shown in following code:
 *
 *		W = rand(m,k);
 *		H = rand(k,n);
 *		for iter=1:maxiter
 *			H = H - epsH * gradH;
 *			W = W - epsW * gradW;
 *		end
 *
 *		This routine computes the necessary stepsize epsW. Using this stepsize the next iteration result
 *		of W is calculated.
 *		
 *		The algorithm is described in Lin (2007) (see references [2]).
 *
 *		This routine uses the BLAS "dgemm" routine for matrix-matrix multiplications
 *
 * Arguments:
 *
 * a		in	pointer to matrix to factorise
 *
 * w0		in/out	pointer to factor matrix w
 *			on exit w0 contains the new solution to w
 *
 * h		in	pointer to factor matrix h
 *
 * grad		in/out	pointer to allocated memory for the gradient
 *			on exit grad contains the gradient
 *
 * wn		in	pointer to allocated memory for possible next w
 *
 * aht		in	pointer to allocated memory for A*H'
 *
 * hht		in	pointer to allocated memory for H*H'
 *
 * d		in	pointer to allocated memory for D = A - W*H
 *
 * tempw	in	pointer to allocated memory for a temporary variable
 *
 * m		in, 	first dimension of a and w/w0
 *
 * n		in, 	second dimension of a and h/h0
 *
 * k		in, 	second dimension of w/w0 and first dimension of h/h0
 *
 * tol		in	tolerance value for stopping condition
 *
 * maxiter	in, 	maximal number of iterations to run
 *
 * function value	out, 	number of iterations actually run
 */

int pg_subprob_w(double *a, double *w0, double *h, double *grad, double *wn, double *aht, double *hht, double * d, double * tempw, int m, int n, int k, double tol, int maxiter);



#endif
