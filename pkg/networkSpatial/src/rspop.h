/*
######################################################################
#
# rspop.h
#
# Written by Carter T. Butts <buttsc@uci.edu>
# Last Modified 05/22/09
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/networkSpatial package.
#
# This file contains headers for rspop.c.
#
######################################################################
*/
#ifndef RSPOP_H
#define RSPOP_H


/*INTERNAL FUNCTIONS--------------------------------------------------------*/



/*R-CALLABLE FUNCTIONS------------------------------------------------------*/

void rqhalton_R(double *coord, int *pn, int *pd, double *base, int *start, double *cmin, double *cmax);

#endif
