/*
######################################################################
#
# utils.c
#
# Written by Carter T. Butts <buttsc@uci.edu>
# Last Modified 12/04/10
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/networkSpatial package.
#
# This file contains code for utilities related to the package.
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
#include "utils.h"


/*INTERNAL FUNCTIONS--------------------------------------------------------*/


double bbAreaCalc(double *coord, int dim, int areamode)
/*Calculate the area of a bounding box performing calculations on the geosphere
if needed.  (The ellipsoid is not currently implemented, but some day it might
be.)  If called with areamode=ACSPHERE, then any dimensions beyond the first
two are treated as ordinary rectilinear extensions.  Coords here should be in
(x0,y0,...,x1,y1,...) form, with the lower coordinates given first.  In the
spherical case, lat/lon order is assumed (and units are taken to be degrees).*/
{
  double area=1.0;
  int i;

  /*If no dimensions, return 0*/
  if(dim<1){
    warning("No dimensions in bbAreaCalc.  Did you intend this?  (Returning 0.0.)");
    return 0.0;
  }

  /*Else, solve for the appropriate area*/
  switch(areamode){
    case ACSPHERE:            /*Rectangular angular region on WGS84 geosphere*/
      if(dim==1){                 /*Assume this is a great circle angle....*/
        area*=geoDistWGS84(0.0,coord[1],0.0,coord[0],GDSPHERE);
      }else{
        area*=(PI/180.0)*WGS84MAJAX*WGS84MAJAX* fabs(sin(coord[dim]*PI/180.0)-sin(coord[0]*PI/180.0))*fabs(coord[dim+1]-coord[1]);
        for(i=2;i<dim;i++)                 /*Treat extra dims as Euclidean*/
          area*=fabs(coord[i+dim]-coord[i]); 
      }
      break;
    case ACEUCLID:            /*Arbitrary rectilinear box*/
      for(i=0;i<dim;i++)
        area*=fabs(coord[i+dim]-coord[i]);
      break;
    default:             /*Shouldn't get here*/
      warning("Illegal area call in bbAreaCalc.  Returning area 0.0.");
      area=0.0;
      break;
  }

  return area;
}


double dSurfaceAreaWGS84(double lat, double lon)
/*Return the surface area element (Jacobian) for the WGS84 ellipsoid at the
specified lat/lon coordinates.  The formula used is

        a^2 (1-e^2) cos phi
  dA = --------------------- dphi dlambda
        1 - e^2 sin^2 phi

where a is the semimajor axis, and e is the eccentricity.
*/
{
  double num,denom,phi,sphi;
  
  phi=lat/180.0*PI;
  num=WGS84MAJAX*WGS84MAJAX*WGS841MESQ*cos(phi);
  sphi=sin(phi);
  denom=1.0-WGS84ESQ*sphi*sphi;
  
  return num/denom;
}


int intIncGeo(int start, double p)
/*
Return start+rgeom(p), with the result forced to be an integer.  If the result
would be too large, it is returned as INT_MAX; this is checked pretty
extensively to ensure against overflow resulting from either the random number
or the addition with the starting point.  This is a special-purpose function,
generally used in random graph generation.
*/
{
  double x;
  
  /*If called with p=0, return the maximum integer.*/
  if(p<=0.0)
    return INT_MAX;
  
  x=rgeom(p);                                 /*Take draw*/
  if(DBL_MAX-x<=(double)start)                /*If too big, return INT_MAX*/
    return INT_MAX;  
  if(x+(double)start>=(double)INT_MAX)
    return INT_MAX;

  return start+(int)x;
}


int isInPoly(double x, double y, double *poly, int n)
/*
Direct adaptation of pnpoly by W. Randolph Franklin, from the comp.graphics.algorithms FAQ.  *poly must be an n x 2 coordinate matrix, with (x,y) being the point to be tested.
*/
{
  int i,j,cross;
  
  cross=0;
  for(i=0,j=n-1;i<n;j=i++){
    if( ((poly[i+n]>y) != (poly[j+n]>y)) &&
	 (x < (poly[j]-poly[i]) * (y-poly[i+n]) / (poly[j+n]-poly[i+n]) + poly[i]) )
       cross = !cross;
  }
  return cross;
}


double geoDistWGS84(double lat0,double lon0,double lat1,double lon1, int mode)
/*For mode=GDELLIPSE, calculate the (geodesic) arc distance between (lat0,lon0) and (lat1,lon1) on the ellipsoid with properties in WGS84MAJAX, WGS84MINAX, and WGS84FLAT (i.e. the WGS-84 Earth ellipsoid).  The algorithm used is based on Vincenty (1975).

Check: calling with (37.95103, 144.4249), (37.65282,143.9265) should return 54974.063.  More examples:
  (46.45938, -109.1177), (37.65282, 143.9265) -> 8195129.041
  (46.45938, -109.1177), (40.0, -107.0) -> 737873.701

Alternately, for mode set to GDSPHERE, a reference sphere is used instead.  This is faster, but sacrifices accuracy.

*/
{
  double rlat0,rlat1,dlon,dlon2,cosrl0,cosrl1,sinrl0,sinrl1,cosdl,sindl;
  double sinsig=0.0,cossig=0.0,sigma=0.0,sina,cossqa=0.0,cos2sm=0.0;
  double C,L=(lon1-lon0)/180.0*PI,B,A,dsigma,usq,rlon0=0.0,rlon1=0.0,gd=0.0;
  
  switch(mode){
    case GDELLIPSE:
      /*Compute the "reduced latitudes"*/
      rlat0=atan((1.0-WGS84FLAT)*tan(lat0/180.0*PI));
      rlat1=atan((1.0-WGS84FLAT)*tan(lat1/180.0*PI));

      /*Define the initial longitude difference*/
      dlon=L;
      dlon2=2*PI;
  
      /*Quick check: return 0 if the points are the same*/
      if((dlon==0.0)&&(lat0=lat1))
        return 0.0;
  
      /*Find an approximating auxillary sphere*/
      while(fabs(dlon-dlon2)>1E-12){
        sinrl0=sin(rlat0);                 /*Precompute angular properties*/
        sinrl1=sin(rlat1);
        sindl=sin(dlon);
        cosrl0=cos(rlat0);
        cosrl1=cos(rlat1);
        cosdl=cos(dlon);
        /*Find approximate angular distance between input points*/
        sinsig=sqrt((cosrl1*sindl)*(cosrl1*sindl) + (cosrl0*sinrl1-sinrl0*cosrl1*cosdl)*(cosrl0*sinrl1-sinrl0*cosrl1*cosdl));
        cossig=sinrl0*sinrl1+cosrl0*cosrl1*cosdl;
        sigma=atan2(sinsig,cossig);
        /*Update the auxillary sphere (and longitude difference thereon)*/
        sina=cosrl0*cosrl1*sindl/sinsig;
        cossqa=1.0-sina*sina;
        if(fabs(cossqa)<1E-12)
          cos2sm=0.0;
        else
          cos2sm=cossig-2*sinrl0*sinrl1/cossqa;
        C=WGS84FLAT/16.0*cossqa*(4.0+WGS84FLAT*(4.0-3.0*cossqa));
        dlon2=dlon;
        dlon=L+(1.0-C)*WGS84FLAT*sina*(sigma+C*sinsig*(cos2sm+ C*cossig*(2.0*cos2sm*cos2sm-1.0)));
      }
      /*Find and save the distance*/
      usq=cossqa*(WGS84MAJAX*WGS84MAJAX-WGS84MINAX*WGS84MINAX) / (WGS84MINAX*WGS84MINAX);
      A=1+usq/16384.0*(4096.0+usq*(usq*(320.0-175.0*usq)-768.0));
      B=usq/1024.0*(256.0+usq*(usq*(74.0-47.0*usq)-128.0)); 
      dsigma=B*sinsig*(cos2sm+B/4.0*(cossig*(2.0*cos2sm*cos2sm-1.0)- B/6.0*cos2sm*(4.0*sinsig*sinsig-3.0)*(4.0*cos2sm*cos2sm-3.0)));
      gd=WGS84MINAX*A*(sigma-dsigma);
      break;
    case GDSPHERE:
      rlat0=lat0/180*PI;
      rlat1=lat1/180*PI;
      rlon0=lon0/180*PI;
      rlon1=lon1/180*PI;
      gd=WGS84MAJAX*acos(cos(rlat0)*cos(rlat1)*cos(rlon0-rlon1)+ sin(rlat0)*sin(rlat1));
      break;
  }

  /*Return the distance*/
  return gd;
}


double maxWGS84GeoDist(double *bb1, double *bb2, int mode)
/*Return the maximum geodesic distance between two lat/lon bounding boxes; this is frankly somewhat ad hoc at present, so it should only be used for crude bounds.*/
{
  double d=0.0;

  /*Should always be an interaction among opposite corners so long as the*/
  /*bounding boxes are not too big (e.g., if they wrap around a pole, this*/
  /*won't work).  For now, I'm ignoring the pathological cases, but this*/
  /*needs to be fixed in the long run.*/
  d=MAX(d,geoDistWGS84(bb1[2],bb1[3],bb2[0],bb2[1],mode));
  d=MAX(d,geoDistWGS84(bb1[0],bb1[1],bb2[2],bb2[3],mode));
  d=MAX(d,geoDistWGS84(bb1[2],bb1[1],bb2[0],bb2[3],mode));
  d=MAX(d,geoDistWGS84(bb1[0],bb1[3],bb2[2],bb2[1],mode));
  
  return d;
}


double minWGS84GeoDist(double *bb1, double *bb2, int mode)
/*Return the minimum geodesic distance between two lat/lon bounding boxes; this is frankly somewhat ad hoc at present, so it should only be used for crude bounds.*/
{
  double d=R_PosInf,lat0=0.0,lat1=0.0,lon0=0.0,lon1=0.0;

  if((bb1[2]>=bb2[0])&&(bb1[0]<=bb2[2])){  /*Latitudes overlap*/
    if((bb1[3]>=bb2[1])&&(bb1[1]<=bb2[3])){  /*Longitudes overlap*/
      /*Boxes intersect - min dist is zero!*/    
      d=0.0;
    }else{                                   /*Longitudes don't overlap*/
      /*Min dist will be const lat line at "top" of boxes for pos lat,*/
      /*or "bottom" for negative lat (have to check for most extreme).*/
      /*(Yes, this is a spherical approximation.)*/
      lat0=MIN(bb1[2],bb2[2]);
      lat1=MAX(bb1[0],bb2[0]);
      if(lat0>fabs(lat1))
        lat1=lat0;
      else
        lat0=lat1;
      d=MIN(geoDistWGS84(lat0,bb1[3],lat1,bb2[1],mode), geoDistWGS84(lat0,bb1[1],lat1,bb2[3],mode));
    }
  }else{                                   /*Latitudes don't overlap*/
    if((bb1[3]>=bb2[1])&&(bb1[1]<=bb2[3])){  /*Longitudes overlap*/
      /*Min dist will be a const long line anywhere (doesn't matter)*/
      lon0=lon1=MIN(bb1[3],bb2[3]);
      d=MIN(geoDistWGS84(bb1[2],lon0,bb2[0],lon1,mode), geoDistWGS84(bb1[0],lon0,bb2[2],lon1,mode));
    }else{                                   /*Longitudes don't overlap*/
      /*Min dist will be between corners; just try 'em all*/
      d=MIN(d,geoDistWGS84(bb1[2],bb1[3],bb2[0],bb2[1],mode));
      d=MIN(d,geoDistWGS84(bb1[0],bb1[1],bb2[2],bb2[3],mode));
      d=MIN(d,geoDistWGS84(bb1[2],bb1[1],bb2[0],bb2[3],mode));
      d=MIN(d,geoDistWGS84(bb1[0],bb1[3],bb2[2],bb2[1],mode));
    }
  }

  return d;
}


double minMinkBBDist(double *bb1, double *bb2, int dim, double me)
/*Return the minimum distance between bounding boxes in a rectilinear space (based on a Minkowski metric).  me here is the Minkowski exponent.  Note that we assume each bb to be of the form (x0,y0,...,x1,y1,...), with length of each series equal to dim.*/
{
  int i;
  double dist=0.0;
  
  /*Work through each dimension*/
  for(i=0;i<dim;i++){
    if(!((bb1[i+dim]>=bb2[i])&&(bb1[i]<=bb2[i+dim]))){ /*Pos dist on dim i?*/
      dist+=pow(MIN(fabs(bb1[i+dim]-bb2[i]),fabs(bb1[i]-bb2[i+dim])),me);
    }
  }

  /*Return the result*/
  return pow(dist,1.0/me);
}


/*List/stack functions*/

element *listInsert(element *head, double val, void *dp)
/*Add a new element to a sorted list, returning a pointer to the updated
list.*/
{
  element *elem,*ep;

  /*Initialize the element*/
  elem=(element *)R_alloc(1,sizeof(struct elementtype));
  elem->val=val;
  elem->dp=dp;
  elem->next=NULL;


  if(head==NULL){  /*If this is the only element, make it the head*/
    return elem;
  }else if(head->val>val){  /*If this is first, make it the head*/
    elem->next=head;
    return elem;
  }else{          /*Otherwise, traverse until we get to the right spot*/
    for(ep=head;(ep->next!=NULL)&&(ep->next->val<val);ep=ep->next);
    if(ep->next==NULL){   /*We ran out of list, apparently*/
      ep->next=elem;
      return head;
    }else{                /*We need to add elem after ep*/
      elem->next=ep->next;
      ep->next=elem;
      return head;
    }
  }
}


element *push(element *head, double val, void *dp)
/*Push a new element to a stack, returning a pointer to the updated
stack.*/
{
  element *elem;

  /*Initialize the element, and place at the top of the stack*/
  elem=(element *)R_alloc(1,sizeof(struct elementtype));
  elem->val=val;
  elem->dp=dp;
  elem->next=head;

  return elem;
}


/*R-Tree functions*/

Rtree *RtreeInsert(Rtree *head, double *bb, int dim, int num, double val, void *dp, int areamode)
/*Insert a new element in the R-Tree pointed to by head, having the specified 
bounding box, value, number, and data pointer.  A pointer to the updated tree 
is  returned.  (To root an empty tree, call with head=NULL.)*/
{
  int i;
  Rtree *new;

  /*Allocate the new element*/  
  new=(Rtree *)R_alloc(1,sizeof(Rtree));
  new->bb=(double *)R_alloc(dim*2,sizeof(double));
  for(i=0;i<2*dim;i++)
    new->bb[i]=bb[i];
  new->dim=dim;
  new->val=val;
  new->num=num;
  new->area=bbAreaCalc(bb,dim,areamode);
  new->dp=dp;
  new->l=NULL;
  new->r=NULL;

  //Rprintf("Inserting new element %d into R-tree\n",num);
  //if(head==NULL)
  //  Rprintf("\tList is empty, creating new root\n");
  
  /*If passed an empy list, just return a pointer to new*/
  if(head==NULL)
    return new;

  //Rprintf("\tList non-empty, recursing....\n");
  /*Otherwise, recursively add, and return a pointer to head*/
  RtreeInsertRecurse(head,new,areamode);
  return head;
}


void RtreeInsertRecurse(Rtree *head, Rtree *new, int areamode)
/*Recursive insertion for R-Trees.  The head node here should be the current
subtree head, and new the current node to be inserted.  The basic cases to be
considered are as follows:

case 0: make current head, new children
  condition: head is a leaf node
  action: create new internal node.  Existing head becomes child of new node.
    new becomes other child.

case 1: merge l,r, add new as child
  condition: bb(new) + current bb is minimal
  action: create new internal node.  Existing head becomes child of new node.
    new becomes other child.

case 2: add new to l tree
  condition: bb(l,new) + bb(r) is minimal
  action: recurse w/head=l, new=new
    
case 3: add new to r tree
  condition: bb(l) + bb(r,new) is minimal
  action: recurse w/head=r, new=new
*/
{
  Rtree *intnode;
  int dim,i;
  double *newbb,amergelr,amergenl,amergenr,mintota;
  
  dim=new->dim;
  
  //Rprintf("\tTrying to place node %d\n",new->num);
  
  /*Case 0: this is a leaf node*/
  if(head->l==NULL){
    //Rprintf("\t\tCase 0: head is a leaf node (%d); creating internal w/leaves %d,%d\n",head->num,new->num,head->num);
    intnode=(Rtree *)R_alloc(1,sizeof(Rtree));     /*Make new internal node*/
    intnode->dim=dim;
    intnode->bb=(double *)R_alloc(2*dim,sizeof(double));
    for(i=0;i<dim;i++){                            /*Copy head into intnode*/
      intnode->bb[i]=head->bb[i];
      intnode->bb[i+dim]=head->bb[i+dim];
    }
    intnode->area=head->area;
    intnode->val=head->val;
    intnode->num=head->num;
    intnode->dp=head->dp;
    intnode->l=head->l;
    intnode->r=head->r;
    head->l=new;                                 /*Make head merge point*/
    head->r=intnode;
    head->dp=NULL;
    head->val=intnode->val+new->val;
    head->num=-1;
    for(i=0;i<dim;i++){
      head->bb[i]=MIN(intnode->bb[i],new->bb[i]);
      head->bb[i+dim]=MAX(intnode->bb[i+dim],new->bb[i+dim]);
    }
    head->area=bbAreaCalc(head->bb,dim,areamode);
  }else{                        /*Not 0 -- depends on bb areas*/
    dim=new->dim;
    if((dim!=head->dim)||(dim!=head->l->dim)||(dim!=head->r->dim))
      error("Incomensurate bounding box dimensions in RtreeInsert.  Exiting.\n");
    newbb=(double *)R_alloc(2*dim,sizeof(double));
    amergelr=head->area+new->area;  /*Easy case: merge l and r*/
    /*Merge new and l, versus r*/
    for(i=0;i<dim;i++){
      newbb[i]=MIN(new->bb[i],head->l->bb[i]);
      newbb[i+dim]=MAX(new->bb[i+dim],head->l->bb[i+dim]);
    }
    amergenl=bbAreaCalc(newbb,dim,areamode);
    /*Merge new and r, versus l*/
    for(i=0;i<dim;i++){
      newbb[i]=MIN(new->bb[i],head->r->bb[i]);
      newbb[i+dim]=MAX(new->bb[i+dim],head->r->bb[i+dim]);
    }
    amergenr=bbAreaCalc(newbb,dim,areamode);
    /*Proceed based on the smallest total area*/
    mintota=MIN(MIN(amergelr,amergenl+head->r->area),amergenr+head->l->area);
    //Rprintf("\tAreas are - (l,r)/n=%f, (n,l)/r=%f, l/(n,r)=%f,min is %f\n",amergelr,amergenl+head->r->area,amergenr+head->l->area,mintota);
    if(amergelr==mintota){                         /*Case 1*/
      //Rprintf("\t\tCase 1: merge current children (%d,%d), set %d as new leaf\n",head->l->num,head->r->num,new->num);
      intnode=(Rtree *)R_alloc(1,sizeof(Rtree));     /*Make new internal node*/
      intnode->dim=dim;
      intnode->bb=newbb;
      for(i=0;i<dim;i++){                            /*Copy head into intnode*/
        intnode->bb[i]=head->bb[i];
        intnode->bb[i+dim]=head->bb[i+dim];
      }
      intnode->area=head->area;
      intnode->val=head->val;
      intnode->num=head->num;
      intnode->dp=head->dp;
      intnode->l=head->l;
      intnode->r=head->r;
      head->l=new;                                 /*Make head merge point*/
      head->r=intnode;
      head->dp=NULL;
      head->val=intnode->val+new->val;
      head->num=-1;
      for(i=0;i<dim;i++){
        head->bb[i]=MIN(intnode->bb[i],new->bb[i]);
        head->bb[i+dim]=MAX(intnode->bb[i+dim],new->bb[i+dim]);
      }
      head->area=bbAreaCalc(head->bb,dim,areamode);
    }else if(amergenl+head->r->area==mintota){     /*Case 2*/
      //Rprintf("\t\tCase 2: add %d to l tree, recurse\n",new->num);
      head->val+=new->val;                           /*Update head*/
      for(i=0;i<dim;i++){
        head->bb[i]=MIN(head->bb[i],new->bb[i]);
        head->bb[i+dim]=MAX(head->bb[i+dim],new->bb[i+dim]);
      }
      head->area=bbAreaCalc(head->bb,dim,areamode);
      RtreeInsertRecurse(head->l,new,areamode);      /*Recurse on l*/
    }else{                                         /*Case 3*/
      //Rprintf("\t\tCase 3: add %d to r tree, recurse\n",new->num);
      head->val+=new->val;                           /*Update head*/
      for(i=0;i<dim;i++){
        head->bb[i]=MIN(head->bb[i],new->bb[i]);
        head->bb[i+dim]=MAX(head->bb[i+dim],new->bb[i+dim]);
      }
      head->area=bbAreaCalc(head->bb,dim,areamode);
      RtreeInsertRecurse(head->r,new,areamode);      /*Recurse on r*/
    }
  }
}


/*R-CALLABLE FUNCTIONS------------------------------------------------------*/

void bbAreaCalc_R(double *area, double *coord, int *dim, int *areamode)
/*Wrapper for bbAreaCalc.  Inputs should be a pointer to the output, and
otherwise should match bbAreaCalc.*/
{
  *area=bbAreaCalc(coord,*dim,*areamode);
}


void geoDistWGS84_R(int *n, double *gd, double *lat0,double *lon0,double *lat1,double *lon1, int *mode)
/*Wrapper for geoDistWGS84.  Inputs should be vectors of length *n, and *mode an integer indicating the computation mode to use (see the .h file).*/
{
  int i;
  
  for(i=0;i<*n;i++){
    gd[i]=geoDistWGS84(lat0[i],lon0[i],lat1[i],lon1[i],*mode);
  }
}

