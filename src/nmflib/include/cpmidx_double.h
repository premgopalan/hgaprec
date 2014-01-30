#ifndef UNIVIE_NMFLIB_CMPIDX_DOUBLE_H
#define UNIVIE_NMFLIB_CMPIDX_DOUBLE_H


/** cpmidx_double - compares two parameters of type idx_double, used for qsort sorting in descending order
 *
 * Purpose:
 *		Defines the comparison function needed by  c standard libraries "qsort()".
 *
 * Arguments:
 *
 * p1			in, 	first parameter to compare
 *
 * p2			in, 	second parameter to compare
 *
 * return value		out, 	returns -1 if p1 < p2, 1 if p1 > p2 and 0 if p1 = p2
 */

int cmpidx_double(const void *p1, const void *p2);


#endif
