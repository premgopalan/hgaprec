#ifndef UNIVIE_NMFLIB_RANDNUMBER_H
#define UNIVIE_NMFLIB_RANDNUMBER_H

/**
 * randnumber - returns a random number between min and max
 *
 * Purpose:
 *		Creating random numbers for initialisation of factor matrices W and H
 *
 * Description:
 *		This routine uses the random number generator rand() from the c standard library
 * 
 * Arguments:
 * 
 * min			in, 	lower interval bound for random numbers
 *
 * max			in, 	upper interval bound for random numbers
 *
 * return value		out, 	generated random number
 */

double randnumber(const int min, const int max);



#endif
