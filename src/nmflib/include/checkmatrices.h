#ifndef UNIVIE_NMFLIB_CHECKMATRICES_H
#define UNIVIE_NMFLIB_CHECKMATRICES_H

/**
 * checkMatrices - checks the passed matrices for non-negativity and compatibility of their dimensions
 * 
 * Purpose:
 *		This routine performs basic checks on the matrices used in nmfDriver to secure the conditions for
 *		the applied algorithms are met
 *
 * Description:
 *		This routine secures, that the matrices passed as arguments are non-negative and have compatible
 *		dimensions.
 *		matrix a - m by n
 *		matrix w - m by k
 *		matrix h - k by n
 *
 * Arguments:
 *
 * a		in, 	pointer to matrix a to check
 *
 * w		in, 	pointer to matrix w to check
 *
 * h		in, 	pointer to matrix h to check
 *
 * m		in, 	first dimension of a and w
 *
 * n		in, 	second dimension of a and h
 *
 * k		in, 	second dimension of w and first dimension of h
 */

int checkMatrices(const double * a, const double * w, const double * h, const int m, \
			 const int n, const int k );


#endif

