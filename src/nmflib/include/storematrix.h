#ifndef UNIVIE_NMFLIB_STOREMATRIX_H
#define UNIVIE_NMFLIB_STOREMATRIX_H


/** storeMatrix
* Purpose:
 * 		Stores a matrix to the specified file
 *
 * Description:
 * 		The file has to be a simple ascii file
 *		The first row has exactly one entry which is the number of rows m
 *		The second row has exactly one entry which is the number of columns n
 *		The next m rows contain the elements of the matrix in row-major-order
 *		Every element is separated by a single space
 *
 * Arguments:
 *
 * fileName	in, 	The path to the file to write
 *
 * m		in, 	first dimension of the stored matrix
 *
 * n		in, 	second dimension of the stored matrix
 *
 * matrix	in, 	m*n matrix to store to the file
 */



void storeMatrix(const char *fileName, const int m, const int n, const double *matrix);


#endif
