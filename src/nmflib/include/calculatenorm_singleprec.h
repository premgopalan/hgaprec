#ifndef UNIVIE_NMFLIB_CALCULATENORM_SINGLEPREC_H
#define UNIVIE_NMFLIB_CALCULATENORM_SINGLEPREC_H


/**
 * calculateNorm - calculates the norm of A-W*H
 *
 * Purpose:
 *		Calculates the frobenius norm of matrix d = A-W*H
 *		The frobenius norm is then divided by sqrt(m*n)
 *
 * Description:
 *		This function uses the Lapack routine "dlange"
 *
 * Arguments:
 *
 * a			in, 	pointer to Matrix a
 *
 * w			in, 	pointer to Matrix w
 *
 * h			in, 	pointer to Matrix h
 *
 * d			in/out, pointer to allocated storage space for matrix d=a-w*h
 *				on exit d = a - w*h
 * m			in, 	first dimension of matrix a and w
 *
 * n			in, 	second dimension of matrix a and h
 *
 * k			in, 	approximation factor, second dimension of matrix w, first dimension of matrix h
 *
 * function value	out, 	frobenius norm of A-W*H / sqrt(m*n)
 */


float calculateNorm_singleprec(float *a, float *w, float *h, float * d, int m, int n, int k);



#endif