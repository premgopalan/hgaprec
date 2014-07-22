#ifndef UNIVIE_NMFLIB_PG_SUBPROB_H_H
#define UNIVIE_NMFLIB_PG_SUBPROB_H_H

/**
 * PG_SUBPROB_H calculates the projected gradient step of factor matrix H used in nmf_pg and nmf_alspg
 * 
 * Purpose:
 *		This routine calculates the projected gradient step of factor matrix H used in nmf_pg and nmf_alspg
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
 *		This routine computes the necessary stepsize epsH. Using this stepsize the next iteration result
 *		of H is calculated.
 *		
 *		The algorithm is described in Lin (2007) (see references [2]).
 *
 *		This routine uses the BLAS "dgemm" routine for matrix-matrix multiplications
 *
 * Arguments:
 *
 * a		in	pointer to matrix to factorise
 *
 * w		in	pointer to factor matrix w
 *
 * h0		in/out	pointer to factor matrix h
 *			on exit h0 contains the new solution to h
 *
 * grad		in/out	pointer to allocated memory for the gradient
 *			on exit grad contains the gradient
 *
 * hn		in	pointer to allocated memory for possible next h
 *
 * wta		in	pointer to allocated memory for W'*A
 *
 * wtw		in	pointer to allocated memory for W'*W
 *
 * d		in	pointer to allocated memory for D = A - W*H
 *
 * temph	in	pointer to allocated memory for a temporary variable
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

int pg_subprob_h(double *a, double *w, double *h0, double *grad, double *hn, double *wta, double *wtw, double * d, double * temph, int m, int n, int k, double tol, int maxiter);




#endif
