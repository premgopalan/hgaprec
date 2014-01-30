#ifndef UNIVIE_NMFLIB_LOADMATRIX_H
#define UNIVIE_NMFLIB_LOADMATRIX_H

/** loadMatrix
 *
 * Purpose:
 * 		Loads a matrix from the specified file and allocates memory according to the matrix dimensions
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
 * fileName	in, 	The path to the file to read
 *
 * m		out, 	first dimension of the loaded matrix
 *
 * n		out, 	second dimension of the loaded matrix
 *
 * matrix	out, 	m*n matrix loaded from the file
 */

void loadMatrix(const char *fileName, int *m, int *n, double **matrix);



#endif

