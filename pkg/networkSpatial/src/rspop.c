/*
######################################################################
#
# rspop.c
#
# Written by Carter T. Butts <buttsc@uci.edu>
# Last Modified 02/01/08
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/networkSpatial package.
#
# This file contains code for simulating spatially-embedded
# populations, and other related things.
#
######################################################################
*/

/*INCLUSIONS----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include "rspop.h"


/*INTERNAL FUNCTIONS--------------------------------------------------------*/

/*R-CALLABLE FUNCTIONS------------------------------------------------------*/

void rqhalton_R(double *coord, int *pn, int *pd, double *base, int *start, double *cmin, double *cmax)
/*
Generate a sample of quasi-random (i.e., low-discrepancy) points using a
d-dimensional Halton sequence.  coord should be an nxd vectorized matrix
(to hold the values), with base being a vector of primes to use as the Halton
bases (one per dimension).  The input start gives the sequence number with which
to start.  Typically, this will be 1, but there are various applications for
which one may want to query a specific element in the Halton sequence (in
which case the starting number should be set>1).
*/
{
  int i,j,n,d;
  double snum,digit,inc; 
  n=*pn;
  d=*pd;

  /*Obtain the Halton numbers*/
  for(i=0;i<n;i++){
    for(j=0;j<d;j++){
      inc=base[j];              /*Get inverse of initial subdivision increment*/
      snum=(double)(i+(*start));          /*Get number in sequence (as double)*/
      coord[i+j*n]=0.0;                               /*Initialize output to 0*/
      while(snum>0.0){         /*Reduce the sequence number to 0 (base change)*/
        digit=fmod(snum,base[j]);  /*Extract the next digit of snum under base*/
        coord[i+j*n]+=digit/inc;            /*Increment output w/current digit*/
        snum=(snum-digit)/base[j];                    /*Reduce sequence number*/
        inc*=base[j];                 /*Increase inverse subdivision increment*/
      }
      coord[i+j*n]*=(cmax[j]-cmin[j]);       /*Translate/scale to output range*/
      coord[i+j*n]+=cmin[j];
    }
  }
}
