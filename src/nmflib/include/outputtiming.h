#ifndef UNIVIE_NMFLIB_OUTPUTTIMING_H
#define UNIVIE_NMFLIB_OUTPUTTIMING_H

#include<sys/time.h>

/** outputTiming - Outputs profile timings
 *
 * Purpose:
 *		Calculate time difference between two stopped times and print the resulting timing
 *
 * Description:
 *		The time is printed in miliseconds
 *
 * Arguments:
 *
 * name		in, 	Text preceeding timing results
 *
 * start	in, 	start of the timing interval
 *
 *end		in, 	end of the timing interval
 */

void outputTiming(const char *name, const struct timeval start, struct timeval end);


#endif
