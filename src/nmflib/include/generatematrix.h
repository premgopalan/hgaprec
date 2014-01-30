#ifndef UNIVIE_NMFLIB_GENERATEMATRIX_H
#define UNIVIE_NMFLIB_GENERATEMATRIX_H

/**
 * generateMatrix - generates an initialization of a matrix
 *
 * Purpose
 *		Allocating memory for a matrix and initialising this matrix using the method specified by "init"
 *
 * Description
 *		In this version only the random initialisation is implemented
 *
 * Arguments:
 *
 * m		in, 	first dimension of matrix
 *
 * n		in, 	second dimension of matrix 
 *
 * k		in,	approximation factor
 *
 * init		in, 	type of initialization
 *
 * min		in, 	lower bound of random numbers
 *
 * max		in, 	upper bound of random numbers
 *
 * matrixW	in/out, pointer to memory storing the matrix (m x k)
 *
 * matrixH	in/out, pointer to memory storing the matrix (k x n)
 * 
 * matrixA	in,	pointer to matrix which should be factorised
 *
 * opts		in,	options for SVD calculation
 */


void generateMatrix(const int m, const int n, const int k, init_t init, const int min, const int max, double **matrixW, double **matrixH, double * matrixA, options_t * opts);



#endif