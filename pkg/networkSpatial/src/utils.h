/*
######################################################################
#
# utils.h
#
# Written by Carter T. Butts <buttsc@uci.edu>
# Last Modified 08/30/09
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/networkSpatial package.
#
# This file contains headers for utils.c.
#
######################################################################
*/
#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

/*MACROS--------------------------------------------------------------------*/

#define MIN(a,b) ((a<b) ? a : b)
#define MAX(a,b) ((a>b) ? a : b)
#define TOINT(a) ((int *)(a))

/*Define Earth ellipsoid properties; by default, we're using WGS-84 in meters*/
#define WGS84MAJAX 6378137.0             /*Length of the major axis, in meters*/
#define WGS84MINAX 6356752.3142          /*Length of the minor axis, in meters*/
#define WGS84FLAT  0.00335281067183099   /*Flattening coefficient*/
#define WGS84IFLAT 298.2572229328696     /*Inverse of the flattening coef*/
#define WGS84E     0.0818191908426       /*Eccentricity*/
#define WGS84ESQ   0.00669437999013      /*Eccentricity squared*/
#define WGS841MESQ 0.993305620010        /*1-(Eccentricity squared)*/

/*Distance computation modes*/
#define GDELLIPSE  0                     /*WGS-84 reference ellipsoid*/
#define GDSPHERE   1                     /*Geosphere*/

/*Area computation modes*/
#define ACSPHERE   0                     /*Geosphere*/
#define ACEUCLID   1                     /*Euclidean*/


/*DATA STRUCTURES-----------------------------------------------------------*/

typedef struct elementtype{
   double val;
   void *dp;
   struct elementtype *next;
} element;


/*(Binary) R-Tree structure.  Each node contains a bounding box (possibly of
zero volume) which must contain the respective bounding boxes of all child 
nodes.  For most applications, each node also contains a value (in val) equal to
the total population counts for its subtree.  A numerical value (num) is also
included to keep track of polygon IDs.  Non-leaf nodes also contain pointers to
left and right sub-elements; leaf nodes instead contain a pointer to arbitrary 
data.

The R-Tree is constructed so as to minimize the bounding box areas at each
level.  Depending upon your coordinate system (e.g., rectilinear geometry
versus spherical), you may need to be careful in how you calculate areas
for a given application.*/
typedef struct Rtreetype{
  double *bb,area,val;
  int dim,num;
  struct Rtreetype *l,*r;
  void *dp;
} Rtree;


/*INTERNAL FUNCTIONS--------------------------------------------------------*/

double bbAreaCalc(double *coord, int dim, int areamode);

double dSurfaceAreaWGS84(double lat, double lon);

int intIncGeo(int start, double p);

int isInPoly(double x, double y, double *poly, int n);

double geoDistWGS84(double lat0,double lon0,double lat1,double lon1, int mode);

double maxWGS84GeoDist(double *bb1, double *bb2, int mode);

double minWGS84GeoDist(double *bb1, double *bb2, int mode);

double minMinkBBDist(double *bb1, double *bb2, int dim, double me);

/*List/stack functions*/

element *listInsert(element *head, double val, void *dp);

element *push(element *head, double val, void *dp);

/*R-Tree functions*/

Rtree *RtreeInsert(Rtree *head, double *bb, int dim, int num, double val, void *dp, int areamode);

void RtreeInsertRecurse(Rtree *head, Rtree *new, int areamode);


/*R-CALLABLE FUNCTIONS------------------------------------------------------*/

void geoDistWGS84_R(int *n, double *gd, double *lat0,double *lon0,double *lat1,double *lon1, int *mode);

void bbAreaCalc_R(double *area, double *coord, int *dim, int *areamode);

#endif
