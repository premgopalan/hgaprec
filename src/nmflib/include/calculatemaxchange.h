#ifndef UNIVIE_NMFLIB_CALCULATEMAXCHANGE_H
#define UNIVIE_NMFLIB_CALCULATEMAXCHANGE_H


/**
 * calculateMaxchange - calculates MaxChange of matrices w or h in the current iteration
 *
 * Purpose:
 *		Calculates the maximum change matrix w or h show from last iteration to the current iteration, which
 *		is used as one of the convergence checks in most implemented nmf algorithms
 *
 * Description:	
 *		Maximum change is calculated as max(max(abs(mat-mat0) / (sqrteps+max(max(abs(mat0))))))
 *		Uses Lapack routine "dlange" to compute the absolute maximum element of a matrix
 *
 * Arguments:
 *
 * mat			in, 	pointer to matrix of the current iteration
 *
 * mat0			in, 	pointer to matrix of the last iteration
 *				on exit mat0 = mat0 - mat
 *
 * m			in, 	first dimension of matrices mat, mat0
 *
 * n			in, 	second dimension of matrices mat, mat0
 *
 * sqrteps		in, 	sqareroot of the machine precision epsilon
 *
 * function value	out, 	calculated maxChange
 */

double calculateMaxchange(double * mat, double * mat0, int m, int n, const double sqrteps);


#endif
