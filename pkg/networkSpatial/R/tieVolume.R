#########################################
#
# tieVolume functions written by Carter T Butts <buttsc@uci.edu>
# wrapper functions written by Zack W Almquist <almquist@uci.edu>
#
# Last Modified 09/28/2011
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/networkSpatial package
#
# This file contains routines for random processes related to 
# spatially embedded networks.  This includes both generation of such
# networks, and simple point processes that may be of utility in
# generating realistic extrapolative simulations from GIS information.
#  tieVolume.core
#  tieVolume.rtree
#
######################################################################




##################
## S3 wrappers for tieVolume
##################

tieVolume <- function(poly, ...){
	UseMethod("tieVolume",poly)
	}

tieVolume.default <-function(poly,...)
{
  stop("Must be a Spatial Polygon")
}

tieVolume.SpatialPolygonsDataFrame<-function(poly,pop,poly2=NULL,pop2=NULL,param,model=c("powlaw","atpowlaw","explaw","loglaw","atanlaw","tflaw","fflaw","cplaw"),tol.pa=c(1e-2,1e-8),tol.mc=c(1e-2,1e-8),trunc=1e-10,gd.mode=c("ellipsoid","sphere"),nq.mode=c("halton","uniform"),maxiter=1e5,log=FALSE,method="core",...){

if(method=="core")	
return(tieVolumecore(poly=poly@polygons,pop=pop,poly2=switch(as.numeric(!is.null(poly2))+1,NULL,poly2@polygons),pop2=pop2,param,model=model,tol.pa=tol.pa,tol.mc=tol.mc,trunc=trunc,gd.mode=gd.mode,nq.mode=nq.mode,maxiter=maxiter,log=log))

if(method=="rtree")
return(tieVolumertree(poly=poly@polygons,pop=pop,param,model=model,tol.pa=tol.pa,tol.mc=tol.mc,trunc=trunc,gd.mode=gd.mode,nq.mode=nq.mode,maxiter=maxiter,log=log))

"Not Valid Method"
}


tieVolume.SpatialPolygons<-function(poly,pop,poly2=NULL,pop2=NULL,param,model=c("powlaw","atpowlaw","explaw","loglaw","atanlaw","tflaw","fflaw","cplaw"),tol.pa=c(1e-2,1e-8),tol.mc=c(1e-2,1e-8),trunc=1e-10,gd.mode=c("ellipsoid","sphere"),nq.mode=c("halton","uniform"),maxiter=1e5,log=FALSE,method="core",...){
	tieVolume.SpatialPolygonsDataFrame(poly=poly,pop=pop,poly2=poly2,pop2=pop2,param,model=model,tol.pa=tol.pa,tol.mc=tol.mc,trunc=trunc,gd.mode=gd.mode,nq.mode=nq.mode,maxiter=maxiter,log=log,method=method)
}






###############
## Base Functions
###############

#Compute the expected tie volumes among a set of polygons
tieVolumertree<-function(poly,pop,param,model=c("powlaw","atpowlaw","explaw","loglaw","atanlaw","tflaw","fflaw","cplaw"),tol.pa=c(1e-2,1e-8),tol.mc=c(1e-2,1e-8),trunc=1e-10,gd.mode=c("ellipsoid","sphere"),nq.mode=c("halton","uniform"),maxiter=1e5,log=FALSE){
  #require(Matrix)
  #Identify the model, and computation modes
  model<-switch(match.arg(model),
    powlaw=1,
    atpowlaw=2,
    explaw=3,
    loglaw=4,
    atanlaw=5,
    tflaw=6,
    fflaw=7,
    cplaw=8
  )
  gd.mode<-switch(match.arg(gd.mode),
    ellipsoid=0,
    sphere=1
  )
  nq.mode<-switch(match.arg(nq.mode),
    halton=0,
    uniform=1
  )
  #Extract polygon information
  n<-length(poly)
  lpoly<-vector(mode="list",n)
  lbb<-vector(mode="list",n)
  for(i in 1:n){
    lpoly[[i]]<-slot(slot(poly[[i]],"Polygons")[[1]],"coords")[,2:1]
    lbb[[i]]<-bbox(slot(poly[[i]],"Polygons")[[1]])[2:1,]
  }
  #Calculate tie volumes
  vol<-.Call("tieVolume_R",n,lpoly,lbb,pop,model,param,tol.pa,tol.mc,trunc, gd.mode,nq.mode,maxiter,log)
  tv<-sparseMatrix(i=vol[[1]],j=vol[[2]],x=vol[[3]],dims=c(n,n))
  #Return the results
  tv
}

tieVolumecore<-function(poly,pop,poly2=NULL,pop2=NULL,param,model=c("powlaw","atpowlaw","explaw","loglaw","atanlaw","tflaw","fflaw","cplaw"),tol.pa=c(1e-2,1e-8),tol.mc=c(1e-2,1e-8),trunc=1e-10,gd.mode=c("ellipsoid","sphere"),nq.mode=c("halton","uniform"),maxiter=1e5,log=FALSE){
  #require(Matrix)
  #Identify the model, and computation modes
  model<-switch(match.arg(model),
    powlaw=1,
    atpowlaw=2,
    explaw=3,
    loglaw=4,
    atanlaw=5,
    tflaw=6,
    fflaw=7,
    cplaw=8
  )
  gd.mode<-switch(match.arg(gd.mode),
    ellipsoid=0,
    sphere=1
  )
  nq.mode<-switch(match.arg(nq.mode),
    halton=0,
    uniform=1
  )
  #Calculate tie volumes
  if(is.null(poly2)){
    if(!is.list(poly))
      poly<-list(poly)
    n<-length(pop)
    tvl<-list()
    tvc<-0
    for(i in 1:n){
      #Get points and area for (outer) polygon
      pts1<-slot(slot(poly[[i]],"Polygons")[[1]],"coords")[,2:1]
      area1<-slot(slot(poly[[i]],"Polygons")[[1]],"area")
      bb1<-bbox(slot(poly[[i]],"Polygons")[[1]])[2:1,]
      for(j in i:n){
        #Get points and area for (outer) polygon
        pts2<-slot(slot(poly[[j]],"Polygons")[[1]],"coords")[,2:1]
        area2<-slot(slot(poly[[j]],"Polygons")[[1]],"area")
        bb2<-bbox(slot(poly[[j]],"Polygons")[[1]])[2:1,]
        #Obtain tie volume estimate
        tvol<-.C("tieVolume_single_R",vol=as.double(0),as.double(pts1), as.integer(NROW(pts1)),as.double(bb1),as.double(pop[i]),as.double(area1), as.double(pts2),as.integer(NROW(pts2)),as.double(bb2),as.double(pop[j]), as.double(area2),as.integer(model),as.double(param),as.double(tol.pa), as.double(tol.mc),as.double(trunc),as.integer(gd.mode),as.integer(nq.mode),as.integer(maxiter),as.integer(i==j),as.integer(log), PACKAGE="networkSpatial",NAOK=TRUE)$vol
        #If volume is non-zero, push to the list (really a stack)
        if(tvol>0){
          temp<-list(i,j,tvol,tvl)
          tvl<-temp
          tvc<-tvc+1
          if(i!=j){
            temp<-list(j,i,tvol,tvl)
            tvl<-temp
            tvc<-tvc+1
          }
        }
      }
    }
    #Write volume information to a sparse matrix
    if(tvc==0){
      tv<-Matrix(0,nrow=n,ncol=n)
    }else{
      rowvec<-rep(0,tvc)
      colvec<-rep(0,tvc)
      tvvec<-rep(0,tvc)
      for(i in 1:tvc){
        rowvec[i]<-tvl[[1]]
        colvec[i]<-tvl[[2]]
        tvvec[i]<-tvl[[3]]
        tvl<-tvl[[4]]
      }
      tv<-sparseMatrix(i=rowvec,j=colvec,x=tvvec,dims=c(n,n))
    }
  }else{
    if(!is.list(poly))
      poly<-list(poly)
    if(!is.list(poly2))
      poly<-list(poly2)
    n<-length(pop)
    m<-length(pop2)
    tvl<-list()
    tvc<-0
    for(i in 1:n){
      #Get points and area for (outer) polygon
      pts1<-slot(slot(poly[[i]],"Polygons")[[1]],"coords")[,2:1]
      area1<-slot(slot(poly[[i]],"Polygons")[[1]],"area")
      bb1<-bbox(slot(poly[[i]],"Polygons")[[1]])[2:1,]
      for(j in 1:m){
        #Get points and area for (outer) polygon
        pts2<-slot(slot(poly2[[j]],"Polygons")[[1]],"coords")[,2:1]
        area2<-slot(slot(poly2[[j]],"Polygons")[[1]],"area")
        bb2<-bbox(slot(poly2[[j]],"Polygons")[[1]])[2:1,]
        #Obtain tie volume estimate
        tvol<-.C("tieVolume_single_R",vol=as.double(0),as.double(pts1), as.integer(NROW(pts1)),as.double(bb1),as.double(pop[i]),as.double(area1), as.double(pts2),as.integer(NROW(pts2)),as.double(bb2),as.double(pop2[j]), as.double(area2),as.integer(model),as.double(param),as.double(tol.pa), as.double(tol.mc),as.double(trunc),as.integer(gd.mode),as.integer(nq.mode),as.integer(maxiter),as.integer(0),as.integer(log), PACKAGE="networkSpatial",NAOK=TRUE)$vol
        #If volume is non-zero, push to the list (really a stack)
        if(tvol>0){
          temp<-list(i,j,tvol,tvl)
          tvl<-temp
          tvc<-tvc+1
        }
      }
    }
    #Write volume information to a sparse matrix
    if(tvc==0){
      tv<-Matrix(0,nrow=n,ncol=m)
    }else{
      rowvec<-rep(0,tvc)
      colvec<-rep(0,tvc)
      tvvec<-rep(0,tvc)
      for(i in 1:tvc){
        rowvec[i]<-tvl[[1]]
        colvec[i]<-tvl[[2]]
        tvvec[i]<-tvl[[3]]
        tvl<-tvl[[4]]
      }
      tv<-sparseMatrix(i=rowvec,j=colvec,x=tvvec,dims=c(n,m))
    }
  }
  #Return the results
  if(any(dim(tv)==1))
    as(tv,"vector")
  else
    tv
}
