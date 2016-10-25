/*
######################################################################
#
# rnspatial.h
#
# Written by Carter T. Butts <buttsc@uci.edu>
# Last Modified 12/04/10
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/networkSpatial package.
#
# This file contains headers for rnspatial.c.
#
######################################################################
*/
#ifndef RNSPATIAL_H
#define RNSPATIAL_H

#include "utils.h"

/*MACROS--------------------------------------------------------------------*/

#define POWLAW       1
#define ATPOWLAW     2
#define EXPLAW       3
#define LOGLAW       4
#define ATANLAW      5
#define TFLAW        6
#define FFLAW        7
#define CPLAW        8

#define QMHALTON     0
#define QMUNIF       1


/*INTERNAL FUNCTIONS--------------------------------------------------------*/

element *rnspatial_RTRecurse(Rtree *p1, Rtree *p2, element *el, int *ne, int n, int directed, double dtol, int tiemodel, double *par, int hd, double ht, double me, int nx, double *x, double *psi, double *sav);

double sif(double d, int model, double *par, int aslog);

double tieVolume(double *poly1, int n1, double *bb1, double pop1, double area1, double *poly2, int n2, double *bb2, double pop2, double area2, int model, double *par, double *ltol, double *lctol, double ltrunc, int gdmode, int quadmode, int maxiter, int issame, int aslog);

element *tieVolume_RTRecurse(Rtree *p1, Rtree *p2, element *vol, int *nvol, int npoly, int model, double *par, double *ltol, double *lctol, double ltrunc, int gdmode, int quadmode, int maxiter, double ltieprob);


/*R-CALLABLE FUNCTIONS------------------------------------------------------*/

SEXP rnspatial_grid_R(SEXP g, SEXP asnet, SEXP coord, SEXP nd, SEXP model, SEXP param, SEXP minkexp, SEXP hdim, SEXP htol, SEXP covar, SEXP baseparam, SEXP tol, SEXP eptol, SEXP subdivisions, SEXP region);

SEXP rnspatial_rtree_R(SEXP g, SEXP asnet, SEXP coord, SEXP nd, SEXP model, SEXP param, SEXP minkexp, SEXP hdim, SEXP htol, SEXP covar, SEXP baseparam, SEXP tol);

SEXP tieVolume_R(SEXP sn, SEXP poly, SEXP bb, SEXP spop, SEXP smodel, SEXP spar, SEXP stol, SEXP sctol, SEXP strunc, SEXP sgdmode, SEXP squadmode, SEXP smaxiter, SEXP saslog);

void tieVolume_single_R(double *vol, double *poly1, int *n1, double *bb1, double *pop1, double *area1, double *poly2, int *n2, double *bb2, double *pop2, double *area2, int *model, double *par, double *tol, double *ctol, double *trunc, int *gdmode, int *quadmode, int *maxiter, int *issame, int *aslog);

#endif
