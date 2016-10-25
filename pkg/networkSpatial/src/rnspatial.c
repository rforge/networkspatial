/*
######################################################################
#
# rnspatial.c
#
# Written by Carter T. Butts <buttsc@uci.edu>
# Last Modified 12/04/10
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/networkSpatial package.
#
# This file contains the code for simulating spatially-embedded 
# networks.
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
#include "networkapi.h"
#include "utils.h"
#include "rnspatial.h"
#include "rspop.h"


/*INTERNAL FUNCTIONS--------------------------------------------------------*/

element *rnspatial_RTRecurse(Rtree *p1, Rtree *p2, element *el, int *ne, int n, int directed, double dtol, int tiemodel, double *par, int hd, double ht, double me, int nx, double *x, double *psi, double *sav)
/*Recursive, R-tree based routine for drawing a spatial Bernoulli graph from a set of rectilinear coordinates.  Arguments are as follows:
  p1,p2 - pointers to R-tree nodes to compare
  el - pointer to a stack storing the edgelist output (element values contain
    num1+n*num2)
  ne - pointer to an integer containing the current stack length (this is 
    incremented whenever el is enlarged)
  n - number of vertices in the set (doesn't need to match the number of
    R-tree leaf nodes, since it is used only for encoding/decoding)
  tiemodel - the tie model (SIF form) to be used
  par - parameter vector for the SIF
  hd - an optional "hieght" dimension, to be treated differently when computing
    distances; this is used to crudely simulate artificial elevation
  ht - "ground level" proximity tolerance to use when assessing the height
    dimension (if any); points within this radius are treated as if they are
    within a vertical structure, while those not within must travel to the
    hyperplane first
  me - the Minkowski exponent to use when computing distances
  nx - the number of covariates (if any) to be used in edge probability
    computation
  x - an optional covariate matrix for edge probability computation
  psi - an optional parameter vector for edge probability computation
  sav - pointer to a count of "avoided" dyads (for diagnostic purposes)
*/
{
  int k,dim,vi,vj;
  double dist,eprob,lp,nd;
  
  dim=p1->dim;
  
  R_CheckUserInterrupt();  /*Might need to exit, if taking too long....*/
  if(p1->l==NULL){    /*p1 is a leaf node*/
    if(p2->l==NULL){     /*p2 is a leaf node*/
      /*Two leaf nodes: try to draw the edge(s)*/
      vi=p1->num;
      vj=p2->num;
      if(vi>vj){  /*Only need to visit each dyad once*/
        /*Compute the inter-point distance (via Minkowski metric)*/
        dist=0.0;
        for(k=0;k<dim;k++)
          if(k!=hd)
            dist+=pow(fabs(p1->bb[k]-p2->bb[k]),me); /*Can use bb for coords!*/
        dist=pow(dist,1.0/me);
        if(hd>-1){        /*If height exists, use it*/
          if(dist<ht)        /*If w/in tol, travel "within" structure*/
            dist+=fabs(p1->bb[hd]-p2->bb[hd]);
          else               /*If not w/in tol, travel "down, over, up"*/
            dist+=fabs(p1->bb[hd])+fabs(p2->bb[hd]);
        }
        /*Compute the base edge probability for this dyad, if x given*/
        if(nx>0){
          for(lp=0.0,k=0;k<nx;k++)
            lp+=psi[k]*abs(x[vi+n*k]-x[vj+n*k]);
          eprob=1.0/(1.0+exp(-lp));
        }else
          eprob=1.0;
        /*Now add the spatial effect*/
        eprob*=sif(dist,tiemodel,par,0);
        /*Finally, draw the edge(s)*/
        if(runif(0.0,1.0)<eprob){
          el=push(el,(double)vi+(double)n*vj,NULL);
          (*ne)++;
        }
        if(directed&&(runif(0.0,1.0)<eprob)){
          el=push(el,(double)vj+(double)n*vi,NULL);
          (*ne)++;
        }
      }
    }else{               /*p2 is an internal node*/
      /*One leaf node: see if we need to recurse on each child of p2*/
      /*Check p1, p2->l*/
      if(p2->l->val==1.0){      /*If both leaves, always recurse*/
        if(p1->num>p2->l->num)  /*Don't bother if we'd only ignore it*/
          el=rnspatial_RTRecurse(p1,p2->l,el,ne,n,directed,dtol,tiemodel,par,hd, ht,me,nx,x,psi,sav);
      }else{
        dist=minMinkBBDist(p1->bb,p2->l->bb,dim,me); /*Get min dist*/
        eprob=log(1.0-sif(dist,tiemodel,par,0));     /*Get pr(non edge)*/
        nd=(p1->val)*(p2->l->val);                   /*Get number of dyads*/
        if(eprob*nd<log(dtol))  /*If pr(no edges)<dtol, then recurse*/
          el=rnspatial_RTRecurse(p1,p2->l,el,ne,n,directed,dtol,tiemodel,par,hd, ht,me,nx,x,psi,sav);
        else
          (*sav)+=nd;                                /*Count saved dyads*/
      }
      /*Check p1, p2->r*/
      if(p2->r->val==1.0){      /*If both leaves, always recurse*/
        if(p1->num>p2->r->num)  /*Don't bother if we'd only ignore it*/
          el=rnspatial_RTRecurse(p1,p2->r,el,ne,n,directed,dtol,tiemodel,par,hd, ht,me,nx,x,psi,sav);
      }else{
        dist=minMinkBBDist(p1->bb,p2->r->bb,dim,me); /*Get min dist*/
        eprob=log(1.0-sif(dist,tiemodel,par,0));     /*Get pr(non edge)*/
        nd=(p1->val)*(p2->r->val);                   /*Get number of dyads*/
        if(eprob*nd<log(dtol))  /*If pr(no edges)<dtol, then recurse*/
          el=rnspatial_RTRecurse(p1,p2->r,el,ne,n,directed,dtol,tiemodel,par,hd, ht,me,nx,x,psi,sav);
        else
          (*sav)+=nd;                                /*Count saved dyads*/
      }
    }
  }else{              /*p1 is an internal node*/
    if(p2->l==NULL){     /*p2 is a leaf node*/
      /*One leaf node: see if we need to recurse on each child of p1*/
      /*Check p1->l, p2*/
      if(p1->l->val==1.0){      /*If both leaves, always recurse*/
        if(p1->l->num>p2->num)  /*Don't bother if we'd only ignore it*/
          el=rnspatial_RTRecurse(p1->l,p2,el,ne,n,directed,dtol,tiemodel,par,hd, ht,me,nx,x,psi,sav);
      }else{
        dist=minMinkBBDist(p1->l->bb,p2->bb,dim,me); /*Get min dist*/
        eprob=log(1.0-sif(dist,tiemodel,par,0));     /*Get pr(non edge)*/
        nd=(p1->l->val)*(p2->val);                   /*Get number of dyads*/
        if(eprob*nd<log(dtol))  /*If pr(no edges)<dtol, then recurse*/
          el=rnspatial_RTRecurse(p1->l,p2,el,ne,n,directed,dtol,tiemodel,par,hd, ht,me,nx,x,psi,sav);
        else
          (*sav)+=nd;                                /*Count saved dyads*/
      }
      /*Check p1->r, p2*/
      if(p1->r->val==1.0){      /*If both leaves, always recurse*/
        if(p1->r->num>p2->num)  /*Don't bother if we'd only ignore it*/
          el=rnspatial_RTRecurse(p1->r,p2,el,ne,n,directed,dtol,tiemodel,par,hd, ht,me,nx,x,psi,sav);
      }else{
        dist=minMinkBBDist(p1->r->bb,p2->bb,dim,me); /*Get min dist*/
        eprob=log(1.0-sif(dist,tiemodel,par,0));     /*Get pr(non edge)*/
        nd=(p1->r->val)*(p2->val);                   /*Get number of dyads*/
        if(eprob*nd<log(dtol))  /*If pr(no edges)<dtol, then recurse*/
          el=rnspatial_RTRecurse(p1->r,p2,el,ne,n,directed,dtol,tiemodel,par,hd, ht,me,nx,x,psi,sav);
        else
          (*sav)+=nd;                                /*Count saved dyads*/
      }
    }else{               /*p2 is an internal node*/
      /*No leaf nodes: see if we need to recurse on each child pair of p1,p2*/
      /*Check p1->l, p2->l*/
      if((p1->l->val==1.0)&&(p2->l->val==1.0)){      /*If both leaves, recurse*/
        if(p1->l->num>p2->l->num)  /*Don't bother if we'd only ignore it*/
          el=rnspatial_RTRecurse(p1->l,p2->l,el,ne,n,directed,dtol,tiemodel,par, hd,ht,me,nx,x,psi,sav);
      }else{
        dist=minMinkBBDist(p1->l->bb,p2->l->bb,dim,me); /*Get min dist*/
        eprob=log(1.0-sif(dist,tiemodel,par,0));        /*Get pr(non edge)*/
        nd=(p1->l->val)*(p2->l->val);                   /*Get number of dyads*/
        if(eprob*nd<log(dtol))  /*If pr(no edges)<dtol, then recurse*/
          el=rnspatial_RTRecurse(p1->l,p2->l,el,ne,n,directed,dtol,tiemodel,par, hd,ht,me,nx,x,psi,sav);
        else
          (*sav)+=nd;                                   /*Count saved dyads*/
      }
      /*Check p1->l, p2->r*/
      if((p1->l->val==1.0)&&(p2->r->val==1.0)){      /*If both leaves, recurse*/
        if(p1->l->num>p2->r->num)  /*Don't bother if we'd only ignore it*/
          el=rnspatial_RTRecurse(p1->l,p2->r,el,ne,n,directed,dtol,tiemodel,par, hd,ht,me,nx,x,psi,sav);
      }else{
        dist=minMinkBBDist(p1->l->bb,p2->r->bb,dim,me); /*Get min dist*/
        eprob=log(1.0-sif(dist,tiemodel,par,0));        /*Get pr(non edge)*/
        nd=(p1->l->val)*(p2->r->val);                   /*Get number of dyads*/
        if(eprob*nd<log(dtol))  /*If pr(no edges)<dtol, then recurse*/
          el=rnspatial_RTRecurse(p1->l,p2->r,el,ne,n,directed,dtol,tiemodel,par, hd,ht,me,nx,x,psi,sav);
        else
          (*sav)+=nd;                                   /*Count saved dyads*/
      }
      /*Check p1->r, p2->l*/
      if((p1->r->val==1.0)&&(p2->l->val==1.0)){      /*If both leaves, recurse*/
        if(p1->r->num>p2->l->num)  /*Don't bother if we'd only ignore it*/
          el=rnspatial_RTRecurse(p1->r,p2->l,el,ne,n,directed,dtol,tiemodel,par, hd,ht,me,nx,x,psi,sav);
      }else{
        dist=minMinkBBDist(p1->r->bb,p2->l->bb,dim,me); /*Get min dist*/
        eprob=log(1.0-sif(dist,tiemodel,par,0));        /*Get pr(non edge)*/
        nd=(p1->r->val)*(p2->l->val);                   /*Get number of dyads*/
        if(eprob*nd<log(dtol))  /*If pr(no edges)<dtol, then recurse*/
          el=rnspatial_RTRecurse(p1->r,p2->l,el,ne,n,directed,dtol,tiemodel,par, hd,ht,me,nx,x,psi,sav);
        else
          (*sav)+=nd;                                   /*Count saved dyads*/
      }
      /*Check p1->r, p2->r*/
      if((p1->r->val==1.0)&&(p2->r->val==1.0)){      /*If both leaves, recurse*/
        if(p1->r->num>p2->r->num)  /*Don't bother if we'd only ignore it*/
          el=rnspatial_RTRecurse(p1->r,p2->r,el,ne,n,directed,dtol,tiemodel,par, hd,ht,me,nx,x,psi,sav);
      }else{
        dist=minMinkBBDist(p1->r->bb,p2->r->bb,dim,me); /*Get min dist*/
        eprob=log(1.0-sif(dist,tiemodel,par,0));        /*Get pr(non edge)*/
        nd=(p1->r->val)*(p2->r->val);                   /*Get number of dyads*/
        if(eprob*nd<log(dtol))  /*If pr(no edges)<dtol, then recurse*/
          el=rnspatial_RTRecurse(p1->r,p2->r,el,ne,n,directed,dtol,tiemodel,par, hd,ht,me,nx,x,psi,sav);
        else
          (*sav)+=nd;                                   /*Count saved dyads*/
      }
    }
  }

  /*Finally, return the current edge list pointer*/
  return el;
}


double sif(double d, int model, double *par, int aslog)
/*Compute the spatial interaction function (tie probability by distance) for the specified model, using the vector of parameters in par.*/
{
  double eprob=0.0;
  
  switch(model){             /*Compute SIF by model*/
    case POWLAW:
      if(aslog)
        eprob=log(par[0])-par[2]*log1p(par[1]*d);
      else
        eprob=par[0]/pow(1.0+par[1]*d,par[2]);
      break;
    case ATPOWLAW:
      if(aslog)
        eprob=log(par[0])-logspace_add(0.0,par[2]*log(par[1]*d));
      else
        eprob=par[0]/(1.0+pow(par[1]*d,par[2]));
      break;
    case EXPLAW:
      if(aslog)
        eprob=log(par[0])-par[1]*d;
      else
        eprob=par[0]*exp(-par[1]*d);
      break;
    case LOGLAW:
      if(aslog)
        eprob=log(1.0+par[2])+log(par[0])-logspace_add(0.0, log(par[2])+par[1]*d);
      else
        eprob=(1.0+par[2])*par[0]/(1.0+par[2]*exp(par[1]*d));
      break;
    case ATANLAW:
      if(aslog)
        eprob=log(par[0])+logspace_sub(0.0,log(2.0/PI)+log(atan(par[1]*d)));
      else
        eprob=par[0]*(1.0-2.0/PI*atan(par[1]*d));
      break;
    case TFLAW:
      if(aslog)
        eprob=MAX(MIN(log(par[2]-par[1]*d),log(par[0])),log(par[3]*par[0]));
      else
        eprob=MAX(MIN(par[2]-par[1]*d,par[0]),par[3]*par[0]);
      break;
    case FFLAW:
      if(d<par[1])
        eprob=(aslog>0) ? log(par[0]) : par[0];
      else
        eprob=(aslog>0) ? log(par[2]*par[0]) : par[2]*par[0];
      break;
    case CPLAW:
      eprob=(aslog>0) ? log(par[0]) : par[0];
      break;
    default:
      error("Illegal edge probability model in sif.\n");
  }
  
  return eprob;
}


double tieVolume(double *poly1, int n1, double *bb1, double pop1, double area1, double *poly2, int n2, double *bb2, double pop2, double area2, int model, double *par, double *ltol, double *lctol, double ltrunc, int gdmode, int quadmode, int maxiter, int issame, int aslog)
/*
Quadrature mode is governed by quadmode, geodesic calculation by gdmode.  If using QMUNIF, be sure to call GetRNGstate prior to calling this function; this is not needed for QMHALTON, which is deterministic.  (Well, so is runif, but you know what I mean.)

Bounding box assumed to be xmin,ymin,xmax,ymax form

Point approximation is used when the tie volume ARE can be shown to be less than exp(ltol).

As of right now, no correction is made for unequal sampling due to lat/lon
issues.  This will produce distortion for very large polygons, but is probably
OK for small cases (especially at non-extreme latitudes).
*/
{
  double d,lepmin,lepmax,lv,coord1[2],coord2[2],base0[2],bmin1[2],bmax1[2];
  double bmin2[2],bmax2[2],ltv=R_NegInf,lse=R_NegInf,lep,ld1,ld2,lstv,lstvsq;
  double mx1,my1,mx2,my2,maxvol,base1[2],lepdiff;
  int n[1],dim[1],sn1[1],sn2[1],iter;
  
  /*If either area is empty, just return 0*/
//  Rprintf("Empty area check\n");
  if((pop1==0.0)||(pop2==0.0))
    return 0.0;

  /*Next, check to see if the point-mass approximation can be used, using*/
  /*bounding boxes*/
//  Rprintf("Point approx check\n");
  if(issame)
    maxvol=log(pop1)+log(pop2-1.0)-log(2.0);
  else
    maxvol=log(pop1)+log(pop2);
  if(!issame){
    lepmin=sif(maxWGS84GeoDist(bb1,bb2,GDSPHERE),model,par,1);
    lepmax=sif(minWGS84GeoDist(bb1,bb2,GDSPHERE),model,par,1);
    lepdiff=logspace_sub(lepmax,lepmin);
    if((lepdiff-lepmin<ltol[0])||(lepdiff<ltol[1])){  /*If ARE/AE low, use apx*/
      mx1=(bb1[2]+bb1[0])/2.0;
      my1=(bb1[3]+bb1[1])/2.0;
      mx2=(bb2[2]+bb2[0])/2.0;
      my2=(bb2[3]+bb2[1])/2.0;
      d=geoDistWGS84(mx1,my1,mx2,my2,gdmode);
      ltv=maxvol+sif(d,model,par,1);
      if(ltv>maxvol)
        error("ltv was %f, maxvol was %f (pop=%f,%f); d=%f, sif=%f\n",ltv, maxvol,pop1,pop2,d,sif(d,model,par,1));
      if(ltv<ltrunc)
        ltv=R_NegInf;
      if(aslog)
        return ltv;
      else
        return exp(ltv);
    }
  }    

  /*Hmm, guess not.  In that case, perform a Monte Carlo quadrature, using*/
  /*Halton sequences or uniform random draws (depending on quadmode).*/
//  Rprintf("MC quadrature\n");
  n[0]=1;                          /*Set up the fixed parameters*/
  base0[0]=2.0;
  base0[1]=3.0;
  base1[0]=5.0;
  base1[1]=7.0;
  dim[0]=2;
  bmin1[0]=bb1[0];
  bmin1[1]=bb1[1];
  bmax1[0]=bb1[2];
  bmax1[1]=bb1[3];
  bmin2[0]=bb2[0];
  bmin2[1]=bb2[1];
  bmax2[0]=bb2[2];
  bmax2[1]=bb2[3];
  ld1=log(pop1)-log(area1);
  ld2=log(pop2)-log(area2);
  sn1[0]=1;                        /*Initialize various things*/
  sn2[0]=1;
  iter=1;
  lstv=R_NegInf;
  lstvsq=R_NegInf;
  lse=R_PosInf;
  ltv=0.0;
  /*Iterate until we hit maxiter, relative uncertainty low enough, or the*/
  /*upper bound of our estimate is below our "zero" threshold.*/
  while((iter<30)||((iter<maxiter)&&(lse-ltv>lctol[0])&& (logspace_add(ltv,lse+0.693)>lctol[1]))){
     R_CheckUserInterrupt();
    /*Draw legal coordinates from each polygon*/
    do{
      if(quadmode==QMHALTON){
        rqhalton_R(coord1,n,dim,base0,sn1,bmin1,bmax1);
        sn1[0]++;
      }else{
        coord1[0]=runif(bmin1[0],bmax1[0]);
        coord1[1]=runif(bmin1[1],bmax1[1]);
      }
    }while(!isInPoly(coord1[0],coord1[1],poly1,n1));
    do{
      if(quadmode==QMHALTON){
        rqhalton_R(coord2,n,dim,base1,sn2,bmin2,bmax2);
        sn2[0]++;
      }else{
        coord2[0]=runif(bmin2[0],bmax2[0]);
        coord2[1]=runif(bmin2[1],bmax2[1]);
      }
    }while(!isInPoly(coord2[0],coord2[1],poly2,n2));
    /*Calculate the instantaneous log volume element*/
    d=geoDistWGS84(coord1[0],coord1[1],coord2[0],coord2[1],gdmode);
//    Rprintf("\tCoord: (%f,%f),(%f,%f) -> %f\n",coord1[0],coord1[1],coord2[0], coord2[1],d);
    lep=sif(d,model,par,1);
    lv=lep;   /*Fix later - note, pulled out density term (const. den apx)*/
    /*Accumulate sum and sum of squares*/
    lstv=logspace_add(lstv,lv);
    lstvsq=logspace_add(lstvsq,2.0*lv);
    /*Update mean and standard error estimates*/
    ltv=maxvol+lstv-log(iter);
    if(ltv>maxvol)
      error("ltv was %f, maxvol was %f (pop=%f,%f); d=%f, sif=%f\n",ltv, maxvol,pop1,pop2,d,sif(d,model,par,1));
    lse=logspace_sub(2.0*maxvol+lstvsq-log(iter),2.0*ltv)/2.0- log(iter)/2.0;
//    Rprintf("\tIter %d: ETV: %f(+/-%f); lograt=%f; lsums %f, %fl lpop %f,%f\n",iter,exp(ltv),1.96*exp(lse),lse-ltv, lstv, lstvsq,log(pop1),log(pop2));
    iter++;
  }
  /*Return the volume estimate*/
//  Rprintf("Finished\n");
  if(ltv<ltrunc)
    ltv=R_NegInf;
  if(aslog)
    return ltv;
  else
    return exp(ltv);
}


element *tieVolume_RTRecurse(Rtree *p1, Rtree *p2, element *vol, int *nvol, int npoly, int model, double *par, double *ltol, double *lctol, double ltrunc, int gdmode, int quadmode, int maxiter, double ltieprob)
/*Recursive, R-tree based routine for calculating interregional tie volumes on a WGS84 polygon set.  Arguments are as follows:
  p1,p2 - pointers to R-tree nodes to compare
  vol - pointer to a stack storing the volume output (element values contain
    num1+npoly*num2, dp points to a double with the log expected tie volume)
  nvol - pointer to an integer containing the current stack length (this is 
    incremented whenever vol is enlarged)
  npoly - number of polygons in the set (doesn't need to match the number of
    R-tree leaf nodes, since it is used only for encoding/decoding)
  model - SIF model code
  par - SIF parameters
  ltol,lctol - Flatness/convergence tolerances; see tieVolume
  ltrunc - Tie volumes below this quantity are truncated to zero, and the
    recursion tree is pruned accordingly; appropriate choice of ltrunc is the
    primary way of gaining efficiency through this algorithm
  gdmode - High-quality distance evaluation mode (note the GDSPHERE is still
    used for crude bounds, but gdmode is used otherwise)
  quadmode,maxiter - Quadrature parameters for tieVolume
  ltieprob - If <=0, it is assumed that we are in a "flat" area and using a
    fixed-probability approximation with tie probability equal to exp(ltieprob).
    Otherwise, it is assumed that we are not in a "flat" area (although we
    may identify such, and recursion subtrees will behave accordingly).
At the root, this function should be called with p1=p2=R-tree head, an empty
(NULL) vol, nvol[0]=0, and ltieprob>0.0.  The root call will return a pointer
to the updated volume stack, and will as a side effect have modified nvol to
be equal to the stack length (a time-saver for subsequent memory allocation).
For information on the tie volume calculation itself, try tieVolume; note that
the same point-mass-approximation used there is used here, and under the
same conditions.
*/
{
  double *dp,ltv=0.0,mind,maxd,ltp,lepmin,lepmax,lepdiff,mx1,my1,mx2,my2;

  R_CheckUserInterrupt();  /*Might need to exit, if taking too long....*/
  //Rprintf("\tRecursion: (%0.2f %0.2f %0.2f %0.2f) [%0.0f] vs (%0.2f %0.2f %0.2f %0.2f) [%0.0f]\n",p1->bb[0],p1->bb[1],p1->bb[2],p1->bb[3],p1->val, p2->bb[0],p2->bb[1],p2->bb[2],p2->bb[3],p2->val);
    
  /*First, if called with a valid tieprob, this means that we should*/ 
  /*either go ahead and compute volumes using the point-mass approximation,*/
  /*or else immediately recurse downward.*/
  if(ltieprob<=0.0){
    //Rprintf("\tUsing point-mass approximation (%f)\n",ltieprob);
    if(p1->l==NULL){    /*p1 is a leaf node*/
      if(p2->l==NULL){     /*p2 is a leaf node*/
        //Rprintf("\t\tp1,p2 both leaf nodes\n");
        /*Compute the tie volume using the point-mass approximation*/
        //if(p1->num>=p2->num){
          ltv=ltieprob+log(p1->val)+log(p2->val);
          if(p1->num==p2->num)
            ltv-=log(2.0);
          if(ltv>ltrunc){    /*If large enough, push to the output stack*/
            dp=(double *)R_alloc(1,sizeof(double));
            dp[0]=ltv;
            vol=push(vol,(double)(p1->num)+(double)npoly*(p2->num),(void *)dp);
            nvol[0]++;
          }
        //}
      }else{              /*p2 is an internal node*/
        //Rprintf("\t\tp1 leaf node, p2 internal\n");
        /*If E(tv)>trunc, recurse on p2's children*/
        if(ltieprob+log(p1->val)+log(p2->l->val)>ltrunc)
          vol=tieVolume_RTRecurse(p1,p2->l,vol,nvol,npoly,model,par,ltol,lctol, ltrunc,gdmode,quadmode,maxiter,ltieprob);
        if(ltieprob+log(p1->val)+log(p2->r->val)>ltrunc)
          vol=tieVolume_RTRecurse(p1,p2->r,vol,nvol,npoly,model,par,ltol,lctol, ltrunc,gdmode,quadmode,maxiter,ltieprob);
      }
    }else{              /*p1 is an internal node*/
      if(p2->l==NULL){     /*p2 is a leaf node*/
        //Rprintf("\t\tp1 internal node, p2 leaf\n");
        /*If E(tv)>trunc, recurse on p1's children*/
        if(ltieprob+log(p1->l->val)+log(p2->val)>ltrunc)
          vol=tieVolume_RTRecurse(p1->l,p2,vol,nvol,npoly,model,par,ltol,lctol, ltrunc,gdmode,quadmode,maxiter,ltieprob);
        if(ltieprob+log(p1->r->val)+log(p2->val)>ltrunc)
          vol=tieVolume_RTRecurse(p1->r,p2,vol,nvol,npoly,model,par,ltol,lctol, ltrunc,gdmode,quadmode,maxiter,ltieprob);
      }else{              /*p2 is an internal node*/
        //Rprintf("\t\tp1,p2 both internal nodes\n");
        /*If E(tv)>trunc, recurse on both child sets*/
        if(ltieprob+log(p1->l->val)+log(p2->l->val)>ltrunc)
          vol=tieVolume_RTRecurse(p1->l,p2->l,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
        if(ltieprob+log(p1->l->val)+log(p2->r->val)>ltrunc)
          vol=tieVolume_RTRecurse(p1->l,p2->r,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
        if(ltieprob+log(p1->r->val)+log(p2->l->val)>ltrunc)
          vol=tieVolume_RTRecurse(p1->r,p2->l,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
        if(ltieprob+log(p1->r->val)+log(p2->r->val)>ltrunc)
          vol=tieVolume_RTRecurse(p1->r,p2->r,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
      }
    }
  }else{  /*Things aren't fixed -- compute or recurse as needed*/
    //Rprintf("\tNOT using fixed-point approximation\n");
    if(p1->l==NULL){    /*p1 is a leaf node*/
      if(p2->l==NULL){     /*p2 is a leaf node*/
        //Rprintf("\t\tp1,p2 both leaf nodes\n");
        /*Compute the tie volume the hard way*/
        ltv=tieVolume((double *)(p1->dp)+1,(int)(((double *)(p1->dp))[0]),p1->bb,p1->val,p1->area,(double *)(p2->dp)+1,(int)(((double *)(p2->dp))[0]),p2->bb,p2->val,p2->area,model,par,ltol,lctol,ltrunc,gdmode, quadmode,maxiter,(p1->num)==(p2->num),1);
        if(ltv>ltrunc){    /*If large enough, push to the output stack*/
          dp=(double *)R_alloc(1,sizeof(double));
          dp[0]=ltv;
          vol=push(vol,(double)(p1->num)+(double)npoly*(p2->num),(void *)dp);
          nvol[0]++;
        }
      }else{              /*p2 is an internal node*/
        //Rprintf("\t\tp1 leaf node, p2 internal\n");
        /*For each p2 child, develop tie prob bounds.  If flat enough (and*/
        /*above trunc), recurse w/fixed tie prob.  If above trunc, recurse.*/
        /*Start with p1, p2->l*/
        maxd=maxWGS84GeoDist(p1->bb,p2->l->bb,GDSPHERE);
        mind=minWGS84GeoDist(p1->bb,p2->l->bb,GDSPHERE);
        lepmin=sif(maxd,model,par,1);
        lepmax=sif(mind,model,par,1);
        if(lepmax+log(p1->val)+log(p2->l->val)>ltrunc){ /*Enough vol to care?*/
          lepdiff=logspace_sub(lepmax,lepmin);
          if((lepdiff-lepmin<ltol[0])||(lepdiff<ltol[1])){ /*Const tp approx?*/
            mx1=(p1->bb[2]+p1->bb[0])/2.0;
            my1=(p1->bb[3]+p1->bb[1])/2.0;
            mx2=(p2->l->bb[2]+p2->l->bb[0])/2.0;
            my2=(p2->l->bb[3]+p2->l->bb[1])/2.0;
            ltp=geoDistWGS84(mx1,my1,mx2,my2,gdmode);
            vol=tieVolume_RTRecurse(p1,p2->l,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltp);
          }else{
            vol=tieVolume_RTRecurse(p1,p2->l,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
          }
        }
        /*Now p1, p2->r*/
        maxd=maxWGS84GeoDist(p1->bb,p2->r->bb,GDSPHERE);
        mind=minWGS84GeoDist(p1->bb,p2->r->bb,GDSPHERE);
        lepmin=sif(maxd,model,par,1);
        lepmax=sif(mind,model,par,1);
        if(lepmax+log(p1->val)+log(p2->r->val)>ltrunc){ /*Enough vol to care?*/
          lepdiff=logspace_sub(lepmax,lepmin);
          if((lepdiff-lepmin<ltol[0])||(lepdiff<ltol[1])){ /*Const tp approx?*/
            mx1=(p1->bb[2]+p1->bb[0])/2.0;
            my1=(p1->bb[3]+p1->bb[1])/2.0;
            mx2=(p2->r->bb[2]+p2->r->bb[0])/2.0;
            my2=(p2->r->bb[3]+p2->r->bb[1])/2.0;
            ltp=geoDistWGS84(mx1,my1,mx2,my2,gdmode);
            vol=tieVolume_RTRecurse(p1,p2->r,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltp);
          }else{
            vol=tieVolume_RTRecurse(p1,p2->r,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
          }
        }
      }
    }else{              /*p1 is an internal node*/
      if(p2->l==NULL){     /*p2 is a leaf node*/
        //Rprintf("\t\tp1 internal node, p2 leaf\n");
        /*For each p1 child, develop tie prob bounds.  If flat enough (and*/
        /*above trunc), recurse w/fixed tie prob.  If above trunc, recurse.*/
        /*Start with p1->l, p2*/
        maxd=maxWGS84GeoDist(p1->l->bb,p2->bb,GDSPHERE);
        mind=minWGS84GeoDist(p1->l->bb,p2->bb,GDSPHERE);
        lepmin=sif(maxd,model,par,1);
        lepmax=sif(mind,model,par,1);
        if(lepmax+log(p1->l->val)+log(p2->val)>ltrunc){ /*Enough vol to care?*/
          lepdiff=logspace_sub(lepmax,lepmin);
          if((lepdiff-lepmin<ltol[0])||(lepdiff<ltol[1])){ /*Const tp approx?*/
            mx1=(p1->l->bb[2]+p1->l->bb[0])/2.0;
            my1=(p1->l->bb[3]+p1->l->bb[1])/2.0;
            mx2=(p2->bb[2]+p2->bb[0])/2.0;
            my2=(p2->bb[3]+p2->bb[1])/2.0;
            ltp=geoDistWGS84(mx1,my1,mx2,my2,gdmode);
            vol=tieVolume_RTRecurse(p1->l,p2,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltp);
          }else{
            vol=tieVolume_RTRecurse(p1->l,p2,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
          }
        }
        /*Now p1->r, p2*/
        maxd=maxWGS84GeoDist(p1->r->bb,p2->bb,GDSPHERE);
        mind=minWGS84GeoDist(p1->r->bb,p2->bb,GDSPHERE);
        lepmin=sif(maxd,model,par,1);
        lepmax=sif(mind,model,par,1);
        if(lepmax+log(p1->r->val)+log(p2->val)>ltrunc){ /*Enough vol to care?*/
          lepdiff=logspace_sub(lepmax,lepmin);
          if((lepdiff-lepmin<ltol[0])||(lepdiff<ltol[1])){ /*Const tp approx?*/
            mx1=(p1->r->bb[2]+p1->r->bb[0])/2.0;
            my1=(p1->r->bb[3]+p1->r->bb[1])/2.0;
            mx2=(p2->bb[2]+p2->bb[0])/2.0;
            my2=(p2->bb[3]+p2->bb[1])/2.0;
            ltp=geoDistWGS84(mx1,my1,mx2,my2,gdmode);
            vol=tieVolume_RTRecurse(p1->r,p2,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltp);
          }else{
            vol=tieVolume_RTRecurse(p1->r,p2,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
          }
        }
      }else{              /*p2 is an internal node*/
        //Rprintf("\t\tp1,p2 both internal\n");
        /*For each p1,p2 child, develop tie prob bounds.  If flat enough (and*/
        /*above trunc), recurse w/fixed tie prob.  If above trunc, recurse.*/
        /*Start with p1->l, p2->l*/
        //Rprintf("\t\t\tTrying p1->l,p2->l\n");
        maxd=maxWGS84GeoDist(p1->l->bb,p2->l->bb,GDSPHERE);
        mind=minWGS84GeoDist(p1->l->bb,p2->l->bb,GDSPHERE);
        lepmin=sif(maxd,model,par,1);
        lepmax=sif(mind,model,par,1);
        //Rprintf("\t\t\t\tmind=%.3f, maxd=%.3f, lepmin=%.3f, lepmax=%.3f\n",mind,maxd,lepmin,lepmax);
        //Rprintf("\t\t\t\tVol <= %f+%f+%f=%f; ltrunc=%f\n",lepmax,log(p1->l->val), log(p2->l->val),lepmax+log(p1->l->val)+log(p2->l->val),ltrunc);
        if(lepmax+log(p1->l->val)+log(p2->l->val)>ltrunc){ /*Enough vol?*/
          lepdiff=logspace_sub(lepmax,lepmin);
          if((lepdiff-lepmin<ltol[0])||(lepdiff<ltol[1])){ /*Const tp approx?*/
            mx1=(p1->l->bb[2]+p1->l->bb[0])/2.0;
            my1=(p1->l->bb[3]+p1->l->bb[1])/2.0;
            mx2=(p2->l->bb[2]+p2->l->bb[0])/2.0;
            my2=(p2->l->bb[3]+p2->l->bb[1])/2.0;
            ltp=geoDistWGS84(mx1,my1,mx2,my2,gdmode);
            vol=tieVolume_RTRecurse(p1->l,p2->l,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltp);
          }else{
            vol=tieVolume_RTRecurse(p1->l,p2->l,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
          }
        }
        /*Next p1->l, p2->r*/
        //Rprintf("\t\t\tTrying p1->l,p2->r\n");
        maxd=maxWGS84GeoDist(p1->l->bb,p2->r->bb,GDSPHERE);
        mind=minWGS84GeoDist(p1->l->bb,p2->r->bb,GDSPHERE);
        lepmin=sif(maxd,model,par,1);
        lepmax=sif(mind,model,par,1);
        //Rprintf("\t\t\t\tmind=%.3f, maxd=%.3f, lepmin=%.3f, lepmax=%.3f\n",mind,maxd,lepmin,lepmax);
        //Rprintf("\t\t\t\tVol <= %f+%f+%f=%f; ltrunc=%f\n",lepmax,log(p1->l->val), log(p2->r->val),lepmax+log(p1->l->val)+log(p2->r->val),ltrunc);
        if(lepmax+log(p1->l->val)+log(p2->r->val)>ltrunc){ /*Enough vol?*/
          lepdiff=logspace_sub(lepmax,lepmin);
          if((lepdiff-lepmin<ltol[0])||(lepdiff<ltol[1])){ /*Const tp approx?*/
            mx1=(p1->l->bb[2]+p1->l->bb[0])/2.0;
            my1=(p1->l->bb[3]+p1->l->bb[1])/2.0;
            mx2=(p2->r->bb[2]+p2->r->bb[0])/2.0;
            my2=(p2->r->bb[3]+p2->r->bb[1])/2.0;
            ltp=geoDistWGS84(mx1,my1,mx2,my2,gdmode);
            vol=tieVolume_RTRecurse(p1->l,p2->r,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltp);
          }else{
            vol=tieVolume_RTRecurse(p1->l,p2->r,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
          }
        }
        /*Now, p1->r, p2->l*/
        //Rprintf("\t\t\tTrying p1->r,p2->l\n");
        maxd=maxWGS84GeoDist(p1->r->bb,p2->l->bb,GDSPHERE);
        mind=minWGS84GeoDist(p1->r->bb,p2->l->bb,GDSPHERE);
        lepmin=sif(maxd,model,par,1);
        lepmax=sif(mind,model,par,1);
        //Rprintf("\t\t\t\tmind=%.3f, maxd=%.3f, lepmin=%.3f, lepmax=%.3f\n",mind,maxd,lepmin,lepmax);
        //Rprintf("\t\t\t\tVol <= %f+%f+%f=%f; ltrunc=%f\n",lepmax,log(p1->r->val), log(p2->l->val),lepmax+log(p1->r->val)+log(p2->l->val),ltrunc);
        if(lepmax+log(p1->r->val)+log(p2->l->val)>ltrunc){ /*Enough vol?*/
          lepdiff=logspace_sub(lepmax,lepmin);
          if((lepdiff-lepmin<ltol[0])||(lepdiff<ltol[1])){ /*Const tp approx?*/
            mx1=(p1->r->bb[2]+p1->r->bb[0])/2.0;
            my1=(p1->r->bb[3]+p1->r->bb[1])/2.0;
            mx2=(p2->l->bb[2]+p2->l->bb[0])/2.0;
            my2=(p2->l->bb[3]+p2->l->bb[1])/2.0;
            ltp=geoDistWGS84(mx1,my1,mx2,my2,gdmode);
            vol=tieVolume_RTRecurse(p1->r,p2->l,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltp);
          }else{
            vol=tieVolume_RTRecurse(p1->r,p2->l,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
          }
        }
        /*Finally, p1->r, p2->r*/
        //Rprintf("\t\t\tTrying p1->r,p2->r\n");
        maxd=maxWGS84GeoDist(p1->r->bb,p2->r->bb,GDSPHERE);
        mind=minWGS84GeoDist(p1->r->bb,p2->r->bb,GDSPHERE);
        lepmin=sif(maxd,model,par,1);
        lepmax=sif(mind,model,par,1);
        //Rprintf("\t\t\t\tmind=%.3f, maxd=%.3f, lepmin=%.3f, lepmax=%.3f\n",mind,maxd,lepmin,lepmax);
        //Rprintf("\t\t\t\tVol <= %f+%f+%f=%f; ltrunc=%f\n",lepmax,log(p1->r->val), log(p2->r->val),lepmax+log(p1->r->val)+log(p2->r->val),ltrunc);
        if(lepmax+log(p1->r->val)+log(p2->r->val)>ltrunc){ /*Enough vol?*/
          lepdiff=logspace_sub(lepmax,lepmin);
          if((lepdiff-lepmin<ltol[0])||(lepdiff<ltol[1])){ /*Const tp approx?*/
            mx1=(p1->r->bb[2]+p1->r->bb[0])/2.0;
            my1=(p1->r->bb[3]+p1->r->bb[1])/2.0;
            mx2=(p2->r->bb[2]+p2->r->bb[0])/2.0;
            my2=(p2->r->bb[3]+p2->r->bb[1])/2.0;
            ltp=geoDistWGS84(mx1,my1,mx2,my2,gdmode);
            vol=tieVolume_RTRecurse(p1->r,p2->r,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltp);
          }else{
            vol=tieVolume_RTRecurse(p1->r,p2->r,vol,nvol,npoly,model,par,ltol, lctol,ltrunc,gdmode,quadmode,maxiter,ltieprob);
          }
        }
      }
    }
  }
  
  /*Finally, return the current volume stack pointer*/
  return vol;
}


/*R-CALLABLE FUNCTIONS------------------------------------------------------*/

SEXP rnspatial_grid_R(SEXP g, SEXP asnet, SEXP coord, SEXP nd, SEXP model, SEXP param, SEXP minkexp, SEXP hdim, SEXP htol, SEXP covar, SEXP baseparam, SEXP tol, SEXP eprobtol, SEXP subdivisions, SEXP region)
{
  int n,i,j,k,vi,vj,directed,pc=0,ind,imod,nx,m=0,nsubd,ns,*subdcount,cn,hd;
  int nhdind,*subdlkup,ne,**subdvec,nnec,compcount;
  double *co,*par,dist,me,eprob=1.0,dtol,*x,*psi,lp,*sdbb,*bb,*subdinc,sav=0.0;
  double ht,maxdist,noeprob=0.0,mineprob=1.0,maxheight=0.0,eptol;
  SEXP tail,head,edgeCheck,tailList,headList,namesEval,valsEval,attrdim,attrn;
  element *hp=NULL,*ep,*ep2,**subd,*nonemptycells,*epi,*epj;
  time_t start,cur,mark;

  int (*locptr)(SEXP);  //Local test pointer
  
  /*Check the time*/
  start=time(NULL);
  /*Load the network API*/
  /*Rprintf("Registering functions...\n");*/
  netRegisterFunctions();
  /*Set things up*/
  GetRNGstate();
  Rprintf("Performing setup\n");
  PROTECT(asnet=coerceVector(asnet,INTSXP)); pc++;
  Rprintf("Did asnet conversion\n");
  if(INTEGER(asnet)[0]){            /*Are we working with a network object?*/
    Rprintf("Network case\n");
    if(netNetSize==NULL)
      Rprintf("Apparently, netNetSize lookup failed!\n");
    if(netNetSize_ptr==NULL)
      Rprintf("netNetSize_ptr is null, too\n");
    locptr=(int (*)(SEXP))R_FindSymbol("networkSize", "network", NULL);
    if(locptr==NULL)
      Rprintf("Local attempt also failed\n");
    if(R_FindSymbol("ModelInitialize", "ergm", NULL)==NULL)
      Rprintf("ERGM lookup failed\n");
    else
      Rprintf("Oddly, ergm lookup worked\n");
    if(R_FindSymbol("isAdjacent", "network", NULL)==NULL)
      Rprintf("isAdjacent lookup failed\n");
    else
      Rprintf("Oddly, isAdjacent lookup worked\n");
    n=netNetSize(g);
    Rprintf("Got network size: %d\n",n);
    directed=netIsDir(g);
    Rprintf("Got directedness: %d\n",directed);
  }else{                            /*If not, assume info in an int vector*/
    Rprintf("Int case\n");
    PROTECT(g=coerceVector(g,INTSXP)); pc++;
    Rprintf("Integer coercion worked\n");
    n=INTEGER(g)[0];
    Rprintf("Network size: %d\n",n);
    directed=INTEGER(g)[1];
    Rprintf("Directedness: %d\n",directed);
  }
  Rprintf("Performing setup 1\n");
  PROTECT(nd=coerceVector(nd,INTSXP)); pc++;  /*Number of dimensions*/
  ind=INTEGER(nd)[0];
  if(hdim!=R_NilValue){                 /*Special "height" dimension*/
    PROTECT(hdim=coerceVector(hdim,INTSXP)); pc++;
    hd = (INTEGER(hdim)[0]<=ind) ? INTEGER(hdim)[0]-1 : -1;
  }else
    hd=-1;
  nhdind=ind-(hd>=0);                              /*# Non-height dims*/
  PROTECT(htol=coerceVector(htol,REALSXP)); pc++;  /*"Height" tolerance*/
  ht=REAL(htol)[0];
  Rprintf("Performing setup 2\n");
  PROTECT(minkexp=coerceVector(minkexp,REALSXP)); pc++;  /*Minkowski exponent*/
  me=REAL(minkexp)[0];
  Rprintf("Performing setup 3\n");
  PROTECT(model=coerceVector(model,INTSXP)); pc++;  /*Spatial model code*/
  imod=INTEGER(model)[0];
  Rprintf("Performing setup 4\n");
  PROTECT(coord=coerceVector(coord,REALSXP)); pc++;  /*Vertex coords*/
  co=REAL(coord);
  Rprintf("Performing setup 5\n");
  PROTECT(param=coerceVector(param,REALSXP)); pc++;  /*Spatial parameters*/
  par=REAL(param);
  Rprintf("Performing setup 6\n");
  PROTECT(edgeCheck=allocVector(LGLSXP,1)); pc++;  /*Edge check variable*/
  INTEGER(edgeCheck)[0]=0;
  PROTECT(covar=coerceVector(covar,REALSXP)); pc++;  /*Base rate covariates*/
  x=REAL(covar);
  nx=length(covar)/n;
  PROTECT(baseparam=coerceVector(baseparam,REALSXP)); pc++; /*Base rate params*/
  psi=REAL(baseparam);
  PROTECT(tol=coerceVector(tol,REALSXP)); pc++;  /*Tolerance for edge omission*/
  dtol=REAL(tol)[0];
  PROTECT(eprobtol=coerceVector(eprobtol,REALSXP)); pc++; /*Eprob het. tol*/
  eptol=REAL(eprobtol)[0];
  PROTECT(subdivisions=coerceVector(subdivisions,INTSXP)); pc++; /*# subdiv*/
  nsubd=INTEGER(subdivisions)[0];
  PROTECT(region=coerceVector(region,REALSXP)); pc++;  /*Regional bounding box*/
  bb=REAL(region);  /*BB is of form dimension x min/max*/

  /*Compute subdivisions, and allocate vertices; each cell is given an*/
  /*element stack, containing the vertices within.  This sorting process*/
  /*is order N.*/
  Rprintf("Subdividing vertices for efficient computation\n");
  subdlkup=(int *)R_alloc(nhdind,sizeof(int));     /*Non-height lookup table*/
  if(nhdind==ind){
    for(i=0;i<ind;i++)
      subdlkup[i]=i;
  }else{
    for(i=0;i<hd;i++)
      subdlkup[i]=i;
    for(i=hd+1;i<ind;i++)
      subdlkup[i-1]=i;
  }
//  for(i=0;i<nhdind;i++)
//    Rprintf("%d ",subdlkup[i]);
//  Rprintf("\n");
  subdinc=(double *)R_alloc(nhdind,sizeof(double));/*Get increment on each dim*/
  for(i=0;i<nhdind;i++){
    subdinc[i]=(bb[subdlkup[i]+ind]-bb[subdlkup[i]])/(double)nsubd;
  }
  ns=(int)pow(nsubd,nhdind);                       /*Get # cells*/
  sdbb=(double *)R_alloc(ns*nhdind*2,sizeof(double)); /*ns x nhdind x 2*/
//  Rprintf("Calculating bounding boxes (ns=%d)\n",ns);
  for(i=0;i<ns;i++){  /*Calculate subdivision bounding boxes*/
    k=i;
    for(j=0;j<nhdind;j++){  /*High/low values for each dimension*/
      cn=k%nsubd;
      sdbb[i+j*ns]=bb[subdlkup[j]]+cn*subdinc[j]; /*Low*/
      sdbb[i+j*ns+nhdind*ns]=bb[subdlkup[j]]+(cn+1.0)*subdinc[j]; /*High*/
      k=(int)floor(k/nsubd);
    }
  }
//  Rprintf("Initializing subd stacks\n");
  subd=(element **)R_alloc(ns,sizeof(element *));  /*Initalize subd stacks*/
  subdcount=(int *)R_alloc(ns,sizeof(int));
  for(i=0;i<ns;i++){
    subd[i]=NULL;
    subdcount[i]=0;
  }
//  Rprintf("Sorting vertices\n");
  for(i=0;i<n;i++){  /*Sort vertices into their appropriate cells*/
//    Rprintf("\t%d\n",i);
    k=0;
    for(j=0;j<nhdind;j++)
      k+=((int)floor((co[i+subdlkup[j]*n]-bb[subdlkup[j]])/subdinc[j]))* ((int)pow(nsubd,j));
//    Rprintf("\t\tk=%d, ns=%d\n",ns);
//    Rprintf("\t\tcheck=%d\n",subd[k]==NULL);
    subd[k]=push(subd[k],(double)i,NULL);
//    Rprintf("\t\tCompleted push\n");
    subdcount[k]++;
    /*PROTECT(temp=allocVector(INTSXP,1));
    INTEGER(temp)[0]=k;
    g=netSetVertexAttrib(g,"cell",temp,i+1);
    UNPROTECT(1);*/
    if(hd>-1)
      maxheight=MAX(maxheight,co[i+hd*n]);
  }
  nonemptycells=NULL;
  nnec=0;
  subdvec=(int **)R_alloc(ns,sizeof(int *));  /*Lame as it is, create vectors*/
  for(i=0;i<ns;i++)
    if(subdcount[i]>0){
      nnec++;
      nonemptycells=push(nonemptycells,(double)i,NULL);
      subdvec[i]=(int *)R_alloc(subdcount[i],sizeof(int));
      for(ep=subd[i],j=0;ep!=NULL;ep=ep->next)
        subdvec[i][j++]=(int)ep->val;
    }
  cur=time(NULL);
  Rprintf("\tTotal elapsed time: %0.2f minutes\n",difftime(cur,start)/60);
  /*Draw the graph*/
  mark=time(NULL);            /*Record time*/
  Rprintf("Drawing graph\n");
  compcount=0.0;
  for(epi=nonemptycells;epi!=NULL;epi=epi->next)  /*Move through cells*/
    for(epj=epi;epj!=NULL;epj=epj->next){
      R_CheckUserInterrupt();  /*Allow interrupts, for the impatient*/
      i=(int)(epi->val);
      j=(int)(epj->val);
      if(compcount-10000000.0*floor(compcount/10000000.0)<1.0){ /*ETA*/
        Rprintf("\tEstimated time to task completion: %0.2f minutes\n", (nnec*(nnec-1.0)-compcount)*difftime(time(NULL),mark)/compcount/60.0);
      }
      /*Compute the minimum inter-cell distance (via Minkowski metric)*/
      dist=0.0;
      if(i!=j){
        for(k=0;k<nhdind;k++){
          if(sdbb[i+k*ns]>sdbb[j+k*ns+nhdind*ns]) /*i is above j on dim k*/
            dist+=pow(sdbb[i+k*ns]-sdbb[j+k*ns+nhdind*ns],me);
          else if(sdbb[j+k*ns]>sdbb[i+k*ns+nhdind*ns]) /*j is above i on dim k*/
            dist+=pow(sdbb[j+k*ns]-sdbb[i+k*ns+nhdind*ns],me);
        }
        dist=pow(dist,1.0/me);
      }
/*      Rprintf("\tCluster %d,%d dist=%f\n",i,j,dist);
      for(k=0;k<ind;k++)
        Rprintf("\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n",sdbb[i+k*ns], sdbb[i+k*ns+ind*ns], sdbb[j+k*ns],sdbb[j+k*ns+ind*ns]);*/
      /*Now, bound the probability of having no edges between i and j*/
      eprob=sif(dist,imod,par,0);
      noeprob=pow(1.0-eprob,subdcount[i]*subdcount[j]);
/*      Rprintf("\tPr(no edges)=%f\n",noeprob);*/
      /*Iff the probability of no edges<dtol, then look more closely*/
      if(noeprob<dtol){
        /*Compute the _maximum_ inter-cell distance (via Minkowski metric)*/
        maxdist=0.0;
        for(k=0;k<nhdind;k++){
          if(sdbb[i+k*ns]>sdbb[j+k*ns+nhdind*ns]) /*i is above j on dim k*/
            maxdist+=pow(sdbb[i+k*ns+nhdind*ns]-sdbb[j+k*ns],me);
          else if(sdbb[j+k*ns]>sdbb[i+k*ns+nhdind*ns]) /*j is above i on dim k*/
            maxdist+=pow(sdbb[j+k*ns+nhdind*ns]-sdbb[i+k*ns],me);
          else                                      /*i and j overlap on dim k*/
            maxdist+=pow(MAX(sdbb[i+k*ns+nhdind*ns]-sdbb[j+k*ns], sdbb[j+k*ns+nhdind*ns]-sdbb[i+k*ns]),me);
        }
        maxdist=pow(maxdist,1.0/me);
        if(hd>-1)                       /*If height exists, use worst case*/
          maxdist+=2.0*maxheight;
        /*Now, find the minimum probability of an i,j edge*/
        mineprob=sif(maxdist,imod,par,0);
        /*If eprob-mineprob>=eptol, no choice but to drill down*/
        if(eprob-mineprob>=eptol){  /*Must examine all dyads*/
  //        Rprintf("%f - %f >= %f, having to drill\n",eprob,mineprob,eptol);
          for(ep=subd[i];ep!=NULL;ep=ep->next)    /*Walk the vertices of i...*/
            for(ep2=subd[j];ep2!=NULL;ep2=ep2->next) /*...and of j...*/
              if((i!=j)||((ep->val)>(ep2->val))){  /*...visiting dyads once.*/
                vi=(int)(ep->val);       /*Extract dyad endpoints*/
                vj=(int)(ep2->val);
                /*Compute the inter-point distance (via Minkowski metric)*/
                dist=0.0;
                for(k=0;k<ind;k++)
                  if(k!=hd)
                    dist+=pow(fabs(co[vi+k*n]-co[vj+k*n]),me);
                dist=pow(dist,1.0/me);
                if(hd>-1){        /*If height exists, use it*/
                  if(dist<ht)     /*If w/in tol, travel "within" structure*/
                    dist+=fabs(co[vi+hd*n]-co[vj+hd*n]);
                  else          /*If not w/in tol, travel "down, over, up"*/
                    dist+=fabs(co[vi+hd*n])+fabs(co[vj+hd*n]);
                }
                /*Rprintf("\tDist=%f\n",dist);*/
                /*Compute the base edge probability for this dyad, if x given*/
                if(nx>0){
                  for(lp=0.0,k=0;k<nx;k++)
                    lp+=psi[k]*abs(x[vi+n*k]-x[vj+n*k]);
                  eprob=1.0/(1.0+exp(-lp));
                }else
                  eprob=1.0;
                /*Compute the spatial effect*/
                eprob*=sif(dist,imod,par,0);
                /*Rprintf("\tEdge prob=%f\n",eprob);*/
                /*Draw edges*/
                if(runif(0.0,1.0)<eprob){
                  hp=push(hp,(double)vi+(double)n*vj,NULL);
                  m++;
                }
                if(directed&&(runif(0.0,1.0)<eprob)){
                  hp=push(hp,(double)vj+(double)n*vi,NULL);
                  m++;
                }
              }
        }else{  /*Use constant probability approximation*/
          /*Use the distance between bounding box centers (no height!)*/
          dist=0.0;
          if(i!=j){
            for(k=0;k<nhdind;k++){
              if(sdbb[i+k*ns]>sdbb[j+k*ns]) /*i is above j on dim k*/
                dist+=pow((sdbb[i+k*ns+nhdind*ns]+sdbb[i+k*ns])/2.0- (sdbb[j+k*ns+nhdind*ns]+sdbb[j+k*ns])/2.0,me);
              else 
                dist+=pow((sdbb[j+k*ns+nhdind*ns]+sdbb[j+k*ns])/2.0- (sdbb[i+k*ns+nhdind*ns]+sdbb[i+k*ns])/2.0,me);
            }
            dist=pow(dist,1.0/me);
          }
          /*Get the edge probability - no covariate effects!*/
          eprob=sif(dist,imod,par,0);
          /*Draw the edges*/
          if(i!=j){ /*i and j are different cells.  Start with i->j*/
//          Rprintf("i=%d, j=%d, counts are %d, %d (prod %d)\n",i,j,subdcount[i],subdcount[j],subdcount[i]*subdcount[j]);
//            dist=(double)m;
            for(k=intIncGeo(0,eprob);k<subdcount[i]*subdcount[j];k=intIncGeo(k, eprob)){
              hp=push(hp,(double)subdvec[i][(int)floor(k/subdcount[j])]+ (double)n*subdvec[j][k%subdcount[j]],NULL);
              m++;
              k++;
            }
//          Rprintf("Added %d edges (should have been about %d)\n", m-(int)dist, (int)(eprob*subdcount[i]*subdcount[j]));
            if(directed){   /*If directed, also consider j->i*/
              for(k=intIncGeo(0,eprob);k<subdcount[j]*subdcount[i]; k=intIncGeo(k,eprob)){
                hp=push(hp,(double)subdvec[j][(int)floor(k/subdcount[i])]+ (double)n*subdvec[i][k%subdcount[i]],NULL);
                m++;
                k++;
              }
            }
          }else{ /*Internal ties; doubt this will ever happen, but....*/
            if(directed){  /*Directed case*/
              for(k=intIncGeo(0,eprob);k<subdcount[i]*(subdcount[i]-1); k=intIncGeo(k,eprob)){
                vj=(int)floor(k/(subdcount[i]-1));
                vi=k%(subdcount[i]-1)+(k%(subdcount[i]-1)>(vj-1));
                hp=push(hp,(double)subdvec[i][vi]+ (double)n*subdvec[i][vj],NULL);
                m++;
                k++;
              }
            }else{       /*Undirected case*/
              for(k=intIncGeo(0,eprob);k< (int)(subdcount[i]*(subdcount[i]-1.0)/2.0);k=intIncGeo(k,eprob)){
                vj=subdcount[i]-2-(int)floor(sqrt(0.25+ 2*(choose((double)subdcount[i],2.0)-k-1))-0.5);
                vi=k+1-vj*(subdcount[i]-1)+vj*(vj+1)/2;
                hp=push(hp,(double)subdvec[i][vi]+ (double)n*subdvec[i][vj],NULL);
                m++;
                k++;
              }
            }
          }
        }
      }else
        sav+=((double)subdcount[i])*subdcount[j];
      compcount++;
    }
  Rprintf("Avoided %.0f dyads; fraction avoided was %0.3f.\n",sav, sav/(n*(n-1.0)/2.0));
  cur=time(NULL);
  Rprintf("\tTotal elapsed time: %0.2f minutes\n",difftime(cur,start)/60);

  /*Add edges to the graph*/
  Rprintf("Adding edges to graph\n");
  if(INTEGER(asnet)[0]){   /*Add edges to the network object*/
    PROTECT(tailList=allocVector(VECSXP,m)); pc++;
    PROTECT(headList=allocVector(VECSXP,m)); pc++;
    PROTECT(namesEval=allocVector(VECSXP,m)); pc++;
    PROTECT(valsEval=allocVector(VECSXP,m)); pc++;
    mark=time(NULL);                                /*Record task start time*/
    for(ep=hp,i=0;i<m;i++){
      R_CheckUserInterrupt();                       /*Allow interrupts*/
      if((i+1)%100000==0){                           /*Estimate time remaining*/
        Rprintf("\tEstimated time to task completion: %0.2f minutes\n", ((m-i)*difftime(time(NULL),mark)/i)/60.0);
      }
      PROTECT(tail=allocVector(INTSXP,1));          /*Form head/tail lists*/
      PROTECT(head=allocVector(INTSXP,1));
      INTEGER(tail)[0]=(int)((ep->val)-n*floor((ep->val)/(double)n))+1;
      INTEGER(head)[0]=(int)floor(ep->val/n)+1;
      SET_VECTOR_ELT(tailList,i,tail);              /*Write all into place*/
      SET_VECTOR_ELT(headList,i,head);
      SET_VECTOR_ELT(namesEval,i,R_NilValue);
      SET_VECTOR_ELT(valsEval,i,R_NilValue);
      UNPROTECT(2);
      ep=ep->next;
    }
    g=netAddEdges(g,tailList,headList,namesEval,valsEval,edgeCheck);
  }else{             /*Not a network object - create an sna edgelist*/
    ne=m*(2-directed);    /*Number of edge entries for sna edgelist*/
    PROTECT(g=allocVector(INTSXP,3*ne)); pc++;
    mark=time(NULL);                                /*Record task start time*/
    for(ep=hp,i=0;i<m;i++){                         /*Add the edge entries*/
      R_CheckUserInterrupt();                       /*Allow interrupts*/
      if((i+1)%100000==0){                           /*Estimate time remaining*/
        Rprintf("\tEstimated time to task completion: %0.2f minutes\n", ((m-i)*difftime(time(NULL),mark)/i)/60.0);
      }
      INTEGER(g)[i*(2-directed)]=(int)((ep->val)-n*floor((ep->val)/(double)n))+ 1;
      INTEGER(g)[i*(2-directed)+ne]=(int)floor(ep->val/n)+1;
      INTEGER(g)[i*(2-directed)+2*ne]=1;
      if(!directed){   /*If undirected, write the back-edges*/
        INTEGER(g)[i*2+1]=(int)floor(ep->val/n)+1; 
        INTEGER(g)[i*2+1+ne]=(int)((ep->val)-n*floor((ep->val)/(double)n))+1;
        INTEGER(g)[i*2+1+2*ne]=1;
      }
      ep=ep->next;
    }
    /*Create the dimension attribute (to make this a matrix)*/
    PROTECT(attrdim = allocVector(INTSXP, 2)); pc++;
    INTEGER(attrdim)[0] = ne; 
    INTEGER(attrdim)[1] = 3;
    setAttrib(g, R_DimSymbol, attrdim);
    /*Add the graph size attribute*/
    PROTECT(attrn = allocVector(INTSXP, 1)); pc++;
    INTEGER(attrn)[0] = n;
    setAttrib(g, install("n"), attrn);
  }
  cur=time(NULL);
  Rprintf("\tTotal elapsed time: %0.2f minutes\n",difftime(cur,start)/60);

  /*Finish up and return*/
  PutRNGstate();
  UNPROTECT(pc);
  return g;
}


SEXP rnspatial_rtree_R(SEXP g, SEXP asnet, SEXP coord, SEXP nd, SEXP model, SEXP param, SEXP minkexp, SEXP hdim, SEXP htol, SEXP covar, SEXP baseparam, SEXP tol)
{
  int n,i,j,directed,pc=0,ind,imod,nx,m=0,hd,nhdind,ne;
  double *co,*par,me,dtol,*x,*psi,*bb,sav=0.0,ht;
  Rtree *rt;
  SEXP tail,head,edgeCheck,tailList,headList,namesEval,valsEval,attrdim,attrn;
  element *hp=NULL,*ep;
  time_t start,cur,mark;

  /*Check the time*/
  start=time(NULL);
  /*Load the network API*/
  /*Rprintf("Registering functions...\n");*/
  netRegisterFunctions();
  /*Set things up*/
  GetRNGstate();
  Rprintf("Performing setup\n");
  PROTECT(asnet=coerceVector(asnet,INTSXP)); pc++;
  if(INTEGER(asnet)[0]){            /*Are we working with a network object?*/
    n=netNetSize(g);
    directed=netIsDir(g);
  }else{                            /*If not, assume info in an int vector*/
    PROTECT(g=coerceVector(g,INTSXP)); pc++;
    n=INTEGER(g)[0];
    directed=INTEGER(g)[1];
  }
  /*Rprintf("Performing setup 1\n");*/
  PROTECT(nd=coerceVector(nd,INTSXP)); pc++;  /*Number of dimensions*/
  ind=INTEGER(nd)[0];
  if(hdim!=R_NilValue){                 /*Special "height" dimension*/
    PROTECT(hdim=coerceVector(hdim,INTSXP)); pc++;
    hd = (INTEGER(hdim)[0]<=ind) ? INTEGER(hdim)[0]-1 : -1;
  }else
    hd=-1;
  nhdind=ind-(hd>=0);                              /*# Non-height dims*/
  PROTECT(htol=coerceVector(htol,REALSXP)); pc++;  /*"Height" tolerance*/
  ht=REAL(htol)[0];
  /*Rprintf("Performing setup 2\n");*/
  PROTECT(minkexp=coerceVector(minkexp,REALSXP)); pc++;  /*Minkowski exponent*/
  me=REAL(minkexp)[0];
  /*Rprintf("Performing setup 3\n");*/
  PROTECT(model=coerceVector(model,INTSXP)); pc++;  /*Spatial model code*/
  imod=INTEGER(model)[0];
  /*Rprintf("Performing setup 4\n");*/
  PROTECT(coord=coerceVector(coord,REALSXP)); pc++;  /*Vertex coords*/
  co=REAL(coord);
  /*Rprintf("Performing setup 5\n");*/
  PROTECT(param=coerceVector(param,REALSXP)); pc++;  /*Spatial parameters*/
  par=REAL(param);
  /*Rprintf("Performing setup 6\n");*/
  PROTECT(edgeCheck=allocVector(LGLSXP,1)); pc++;  /*Edge check variable*/
  INTEGER(edgeCheck)[0]=0;
  PROTECT(covar=coerceVector(covar,REALSXP)); pc++;  /*Base rate covariates*/
  x=REAL(covar);
  nx=length(covar)/n;
  PROTECT(baseparam=coerceVector(baseparam,REALSXP)); pc++; /*Base rate params*/
  psi=REAL(baseparam);
  PROTECT(tol=coerceVector(tol,REALSXP)); pc++;  /*Tolerance for edge omission*/
  dtol=REAL(tol)[0];

  /*Create an R-tree for the point set; this should take around O(n log n)*/
  Rprintf("Building spatial R-tree\n");
  rt=NULL;
  for(i=0;i<n;i++){
    bb=(double *)R_alloc(2*ind,sizeof(double));
    for(j=0;j<ind;j++)
      bb[j]=bb[j+ind]=co[i+j*n];
    rt=RtreeInsert(rt,bb,ind,i,1.0,NULL,ACEUCLID);
  }
  cur=time(NULL);
  Rprintf("\tTotal elapsed time: %0.2f minutes\n",difftime(cur,start)/60);


  /*Draw the graph, using recursive double R-tree descent*/
  mark=time(NULL);            /*Record time*/
  Rprintf("Drawing graph\n");
  m=0;
  sav=0.0;
  hp=rnspatial_RTRecurse(rt,rt,NULL,&m,n,directed,dtol,imod,par,hd,ht,me,nx,x, psi,&sav);
  Rprintf("Avoided %.0f ordered dyads; fraction avoided was %0.3f.\n",sav, sav/(n*n));
  cur=time(NULL);
  Rprintf("\tTotal elapsed time: %0.2f minutes\n",difftime(cur,start)/60);

  /*Add edges to the graph*/
  Rprintf("Adding edges to graph\n");
  if(INTEGER(asnet)[0]){   /*Add edges to the network object*/
    PROTECT(tailList=allocVector(VECSXP,m)); pc++;
    PROTECT(headList=allocVector(VECSXP,m)); pc++;
    PROTECT(namesEval=allocVector(VECSXP,m)); pc++;
    PROTECT(valsEval=allocVector(VECSXP,m)); pc++;
    mark=time(NULL);                                /*Record task start time*/
    for(ep=hp,i=0;i<m;i++){
      R_CheckUserInterrupt();                       /*Allow interrupts*/
      if((i+1)%100000==0){                           /*Estimate time remaining*/
        Rprintf("\tEstimated time to task completion: %0.2f minutes\n", ((m-i)*difftime(time(NULL),mark)/i)/60.0);
      }
      PROTECT(tail=allocVector(INTSXP,1));          /*Form head/tail lists*/
      PROTECT(head=allocVector(INTSXP,1));
      INTEGER(tail)[0]=(int)((ep->val)-n*floor((ep->val)/(double)n))+1;
      INTEGER(head)[0]=(int)floor(ep->val/n)+1;
      SET_VECTOR_ELT(tailList,i,tail);              /*Write all into place*/
      SET_VECTOR_ELT(headList,i,head);
      SET_VECTOR_ELT(namesEval,i,R_NilValue);
      SET_VECTOR_ELT(valsEval,i,R_NilValue);
      UNPROTECT(2);
      ep=ep->next;
    }
    g=netAddEdges(g,tailList,headList,namesEval,valsEval,edgeCheck);
  }else{             /*Not a network object - create an sna edgelist*/
    ne=m*(2-directed);    /*Number of edge entries for sna edgelist*/
    PROTECT(g=allocVector(INTSXP,3*ne)); pc++;
    mark=time(NULL);                                /*Record task start time*/
    for(ep=hp,i=0;i<m;i++){                         /*Add the edge entries*/
      R_CheckUserInterrupt();                       /*Allow interrupts*/
      if((i+1)%100000==0){                           /*Estimate time remaining*/
        Rprintf("\tEstimated time to task completion: %0.2f minutes\n", ((m-i)*difftime(time(NULL),mark)/i)/60.0);
      }
      INTEGER(g)[i*(2-directed)]=(int)((ep->val)-n*floor((ep->val)/(double)n))+ 1;
      INTEGER(g)[i*(2-directed)+ne]=(int)floor(ep->val/n)+1;
      INTEGER(g)[i*(2-directed)+2*ne]=1;
      if(!directed){   /*If undirected, write the back-edges*/
        INTEGER(g)[i*2+1]=(int)floor(ep->val/n)+1; 
        INTEGER(g)[i*2+1+ne]=(int)((ep->val)-n*floor((ep->val)/(double)n))+1;
        INTEGER(g)[i*2+1+2*ne]=1;
      }
      ep=ep->next;
    }
    /*Create the dimension attribute (to make this a matrix)*/
    PROTECT(attrdim = allocVector(INTSXP, 2)); pc++;
    INTEGER(attrdim)[0] = ne; 
    INTEGER(attrdim)[1] = 3;
    setAttrib(g, R_DimSymbol, attrdim);
    /*Add the graph size attribute*/
    PROTECT(attrn = allocVector(INTSXP, 1)); pc++;
    INTEGER(attrn)[0] = n;
    setAttrib(g, install("n"), attrn);
  }
  cur=time(NULL);
  Rprintf("\tTotal elapsed time: %0.2f minutes\n",difftime(cur,start)/60);

  /*Finish up and return*/
  PutRNGstate();
  UNPROTECT(pc);
  return g;
}


void tieVolume_single_R(double *vol, double *poly1, int *n1, double *bb1, double *pop1, double *area1, double *poly2, int *n2, double *bb2, double *pop2, double *area2, int *model, double *par, double *tol, double *ctol, double *trunc, int *gdmode, int *quadmode, int *maxiter, int *issame, int *aslog)
/*R wrapper function for tieVolume.  The resulting volume is put in vol; note that tolerances should be passed in non-log form here.*/
{
  double ltol[2],lctol[2];
  
  /*Log the tolerance parameters*/
  ltol[0]=log(tol[0]);
  ltol[1]=log(tol[1]);
  lctol[0]=log(ctol[0]);
  lctol[1]=log(ctol[1]);
  
  GetRNGstate();
  vol[0]=tieVolume(poly1,*n1,bb1,*pop1,*area1,poly2,*n2,bb2,*pop2,*area2, *model,par,ltol,lctol,log(*trunc),*gdmode,*quadmode,*maxiter,*issame,*aslog);
  PutRNGstate();
}


SEXP tieVolume_R(SEXP sn, SEXP poly, SEXP bb, SEXP spop, SEXP smodel, SEXP spar, SEXP stol, SEXP sctol, SEXP strunc, SEXP sgdmode, SEXP squadmode, SEXP smaxiter, SEXP saslog)
/*Compute the tie volumes for all pairs of polygons in poly, with population 
counts in pop.  The result is returned as k x 3 array with row entries poly1,
poly2, and volume.  (This is suitable for conversion to a sparse matrix form.)
Computation is performed by first using an R-Tree to process the polygons, and
then computing volumes only for pairs of polygons whose expected volumes can be
bounded above the truncation threshold.  (Without this, it's simply impossible
to scale these computations to large regions.)

poly - list of polygons
sn - length of poly
bb - list of bounding boxes
spop - vector of population counts
model - code for SIF to use
par - SIF parameters
tol - (unlogged) tolerance parameter vector; first number is fractional
  tolerance for change across a polygon to use the point mass approximation
  (sets maximum diff as a fraction of lower bound), second number is absolute
  error tolerance (i.e., absolute diff that can be ignored, regardless of
  lower bound on total interaction)
ctol - (unlogged) tolerance for convergence; first number is the fractional
  tolerance for quadrature uncertainty (expressed as tolerable standard error
  to mean estimate ratio), second number is a truncation threshold (if upper
  bound on estimate is below this, stop itrerating)
trunc - (unlogged) truncation threshold; estimates generated by any mechanism
  below trunc are returned as 0.0.  This is a useful approximation that greatly
  reduces storage demands (since output is ultimately placed in sparse matrix
  form).
gdmode - geodesic distance computation mode (see .h)
quadmode - quadrature mode
maxiter - maximum iterations to use during quadrature computation
aslog - boolean; return the volume on log scale?
*/
{
  Rtree *rt;
  int i,j,n,areamode,stacklen=0,pc=0,gdmode,quadmode,maxiter,model,aslog;
  double *pop,ltol[2],lctol[2],ltrunc,*par;
  element *volstack,*p;
  SEXP x,y,vol,outl,polytemp,pts;
  
  /*Perform initial coercion and such*/
  //Rprintf("Entering coercion phase\n");
  PROTECT(spop=coerceVector(spop,REALSXP)); pc++;
  pop=REAL(spop);
  //Rprintf("1\n");
  PROTECT(stol=coerceVector(stol,REALSXP)); pc++;
  ltol[0]=log(REAL(stol)[0]);
  ltol[1]=log(REAL(stol)[1]);
  //Rprintf("2\n");
  PROTECT(sctol=coerceVector(sctol,REALSXP)); pc++;
  lctol[0]=log(REAL(sctol)[0]);
  lctol[1]=log(REAL(sctol)[1]);
  //Rprintf("3\n");
  PROTECT(strunc=coerceVector(strunc,REALSXP)); pc++;
  ltrunc=log(REAL(strunc)[0]);
  //Rprintf("4\n");
  PROTECT(spar=coerceVector(spar,REALSXP)); pc++;
  par=REAL(spar);
  //Rprintf("5\n");
  PROTECT(smodel=coerceVector(smodel,INTSXP)); pc++;
  model=INTEGER(smodel)[0];
  PROTECT(smaxiter=coerceVector(smaxiter,INTSXP)); pc++;
  maxiter=INTEGER(smaxiter)[0];
  PROTECT(sgdmode=coerceVector(sgdmode,INTSXP)); pc++;
  gdmode=INTEGER(sgdmode)[0];
  PROTECT(squadmode=coerceVector(squadmode,INTSXP)); pc++;
  quadmode=INTEGER(squadmode)[0];
  areamode=ACSPHERE;
  PROTECT(sn=coerceVector(sn,INTSXP)); pc++;
  n=INTEGER(sn)[0];
  PROTECT(saslog=coerceVector(saslog,INTSXP)); pc++;
  aslog=INTEGER(saslog)[0];
  
  /*Build the R-tree*/
  rt=NULL;
  //Rprintf("Entering R-tree phase\n");
  for(i=0;i<n;i++)
    if(pop[i]>0.0){
      //Rprintf("\tPolygon %d (length=%d)\n",i,length(VECTOR_ELT(poly,i)));
      SET_VECTOR_ELT(bb,i,coerceVector(VECTOR_ELT(bb,i),REALSXP));
      PROTECT(polytemp=coerceVector(VECTOR_ELT(poly,i),REALSXP));
      PROTECT(pts=allocVector(REALSXP,length(polytemp)+1));
      REAL(pts)[0]=length(polytemp)/2;
      for(j=0;j<length(polytemp);j++)
        REAL(pts)[j+1]=REAL(polytemp)[j];
      SET_VECTOR_ELT(poly,i,pts);
      UNPROTECT(2);
      //Rprintf("\t\tInserting...\n");
      rt=RtreeInsert(rt,REAL(VECTOR_ELT(bb,i)),2,i,pop[i], (void *)(REAL(VECTOR_ELT(poly,i))),ACSPHERE);
    }

  /*Now, perform a double-descending traversal of rt to compare polygons*/
  //Rprintf("Entering traversal phase\n");
  GetRNGstate();
  volstack=tieVolume_RTRecurse(rt,rt,NULL,&stacklen,n,model,par,ltol,lctol, ltrunc,gdmode,quadmode,maxiter,2.0);
  PutRNGstate();
  
  /*Now, volstack contains our volume information.  Let's convert it into*/
  /*a more useful vector form (for subsequent processing by Matrix).*/
  //Rprintf("Entering data process phase\n");
  PROTECT(x=allocVector(INTSXP,stacklen)); pc++;
  PROTECT(y=allocVector(INTSXP,stacklen)); pc++;
  PROTECT(vol=allocVector(REALSXP,stacklen)); pc++;
  for(i=0,p=volstack;p!=NULL;i++,p=p->next){
    INTEGER(x)[i]=(int)fmod(p->val,(double)n)+1;
    INTEGER(y)[i]=(int)floor(p->val/n)+1;
    if(aslog)
      REAL(vol)[i]=((double *)(p->dp))[0];
    else
      REAL(vol)[i]=exp(((double *)(p->dp))[0]);
    if(i>=stacklen) //REMOVE THIS LATER.....
      Rprintf("ACK!  stacklen too small in tieVolume!\n");
    else if((i<stacklen-1)&&(p->next==NULL))  //ALSO REMOVE.....
      Rprintf("ACK!  stacklen too large in tieVolume!\n");
  }
  
  //Rprintf("Entering data return phase\n");
  /*Finally, return everything as a list*/
  PROTECT(outl=allocVector(VECSXP,3)); pc++;
  SET_VECTOR_ELT(outl,0,x); 
  SET_VECTOR_ELT(outl,1,y); 
  SET_VECTOR_ELT(outl,2,vol);
  UNPROTECT(pc);
  return outl;
}
