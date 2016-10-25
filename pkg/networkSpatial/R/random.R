######################################################################
#
# random.R
#
# Written by Carter T. Butts <buttsc@uci.edu>.
#
# Last Modified 07/19/10
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/networkSpatial package
#
# This file contains routines for random processes related to 
# spatially embedded networks.  This includes both generation of such
# networks, and simple point processes that may be of utility in
# generating realistic extrapolative simulations from GIS information.
#
# Contents:
#  rnspatial
#  rnspatial.grid
#  rnspatial.rtree
#  rqhalton
#  rspop
#  tieVolume
#  tieVolume.rtree
#
######################################################################


#rnspatial - Draw one or more random network objects, using a spatial model
rnspatial<-function(n,coord,param,accel.method=c("grid","rtree"),...){
  switch(match.arg(accel.method),
    grid=rnspatial.grid(n=n,coord=coord,param=param,...),
    rtree=rnspatial.rtree(n=n,coord=coord,param=param,...)
  )
}


#rnspatial.grid - Draw one or more random network objects, using a spatial model
#(grid implementation)
rnspatial.grid<-function(n,coord,param,model=c("powlaw","atpowlaw","explaw","loglaw","atanlaw","tflaw","fflaw","cplaw"),directed=FALSE,return.as.network=TRUE,minkowski.exp=2,height.dimension=NULL,height.tol=0,covar=NULL,psi=NULL,tol=0.99,epdiff.tol=1e-9,subdivisions=10,...){
  #If drawing multiple graphs, recursively apply
  if(n>1)
    return(lapply(rep(1,n),rnspatial.grid,coord=coord,param=param,model=model, directed=directed,return.as.network=return.as.network, minkowski.exp=minkowski.exp,height.dimension=height.dimension, height.tol=height.tol,covar=covar,psi=psi,tol=tol,subdivisions=subdivisions,...))
  #Find out how many dimensions we have
  if(is.matrix(coord)){
    nd<-NCOL(coord)
    nv<-NROW(coord)
  }else{
    nd<-1
    nv<-length(coord)
  }    
  #Check for psi/covar
  if(!is.null(covar)){
    if(is.null(psi))
      stop("psi must be supplied when covariates are given in rnspatial.grid.")
    else if(length(covar)/nv>length(psi))
      stop("psi length does not match number of covariates in rnspatial.grid.")
  }
  #Find the bounding region for the coordinate set
  if(nd>1){
    region<-t(apply(coord,2,range))
    region[,1]<-region[,1]-1e-8      #Pad the region to avoid numerical issues
    region[,2]<-region[,2]+1e-8
  }else{
    region<-range(coord)
    region<-region+c(-1,1)*1e-8      #Pad the region to avoid numerical issues
  }
  #Begin by initializing an empty graph, or else a vector of characteristics
  if(return.as.network)
    net<-network.initialize(nv,directed=directed)
  else
    net<-c(nv,directed)
  #Identify the model
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
  #Draw and return the graph  
  print(net)
  cat("Entering C space\n")
  .Call("rnspatial_grid_R",net,return.as.network,coord,nd,model,param, minkowski.exp, height.dimension,height.tol,covar,psi,tol,epdiff.tol, subdivisions,region, PACKAGE="networkSpatial")
}


#rnspatial.rtree - Draw one or more random network objects, using a spatial 
#model (R-tree implementation)
rnspatial.rtree<-function(n,coord,param,model=c("powlaw","atpowlaw","explaw","loglaw","atanlaw","tflaw","fflaw","cplaw"),directed=FALSE,return.as.network=TRUE,minkowski.exp=2,height.dimension=NULL,height.tol=0,covar=NULL,psi=NULL,tol=0.99,...){
  #If drawing multiple graphs, recursively apply
  if(n>1)
    return(lapply(rep(1,n),rnspatial.rtree,coord=coord,param=param,model=model, directed=directed,return.as.network=return.as.network, minkowski.exp=minkowski.exp,height.dimension=height.dimension, height.tol=height.tol,covar=covar,psi=psi,tol=tol,...))
  #Find out how many dimensions we have
  if(is.matrix(coord)){
    nd<-NCOL(coord)
    nv<-NROW(coord)
  }else{
    nd<-1
    nv<-length(coord)
  }    
  #Check for psi/covar
  if(!is.null(covar)){
    if(is.null(psi))
      stop("psi must be supplied when covariates are given in rnspatial.rtree.")
    else if(length(covar)/nv>length(psi))
      stop("psi length does not match number of covariates in rnspatial.rtree.")
  }
  #Begin by initializing an empty graph, or else a vector of characteristics
  if(return.as.network)
    net<-network.initialize(nv,directed=directed)
  else
    net<-c(nv,directed)
  #Identify the model
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
  #Draw and return the graph
  .Call("rnspatial_rtree_R",net,return.as.network,coord,nd,model,param, minkowski.exp, height.dimension,height.tol,covar,psi,tol, PACKAGE="networkSpatial")
}


#Not truly pseudorandom, this routine draws quasi-random numbers from a Halton sequence
rqhalton<-function(n,d=1,min=rep(0,d),max=rep(1,d),base=listPrimes(d,start=2),seq.start=1,log=FALSE){
  coord<-.C("rqhalton_R",coord=as.double(matrix(0,n,d)),as.integer(n), as.integer(d), as.double(base),as.integer(seq.start),as.double(min), as.double(max),PACKAGE="networkSpatial")$coord
  if(d>1)
    coord<-matrix(coord,n,d)
  coord
}


#Draw a random spatially embedded population from geographical and household
#data, using either uniform or quasi-random within-block placement and
#stacking.  Note that calculations are in meters here.
rspop<-function(geo,pop,method=c("uniform","halton"),stack.rad=10,stack.dis=4,household.jitter=5,projection.point=NULL){
  require(rgdal)
  #First, fix projection information for raw shape information
  if(is.na(proj4string(geo)))
    proj4string(geo)<-CRS("+proj=longlat +datum=WGS84")
  #Next orthographically project the data about an appropriate point
  if(is.null(projection.point)){
    #Get master bounding box (in long/lat)
    cmin<-c(Inf,Inf)
    cmax<-c(-Inf,-Inf)
    for(i in 1:length(geo[,1])){
      cmin<-pmin(cmin,geo[i,1]@bbox[,1])
      cmax<-pmax(cmin,geo[i,1]@bbox[,2])
    }
    projection.point<-(cmin+cmax)/2           #Get long/lat center of the area
  }  
  #Now, perform the projection
  geoprj<-geo
  geoprj<-spTransform(geoprj,CRS(paste("+proj=ortho +lat_0=",projection.point[2]," +lon_0=",projection.point[1], collapse="",sep="")))
  #Draw points based on the population assignment model
  totpop<-sum(sapply(pop,function(z){sum(z[,1]*z[,2])}))
  x<-rep(NA,totpop)
  y<-rep(NA,totpop)
  z<-rep(NA,totpop)
  pc<-1
  for(i in 1:length(pop)){             #Walk polygons
    cat("Processing polygon",i,"\n")
    shape<-try(as(geoprj[i,1],"owin"))
    if(class(shape)=="try-error"){  #If this failed, deal with sp's shit
      temp<-SpatialPolygons(list(Polygons(list(slot((slot(geoprj[i,1], "polygons")[[1]]),"Polygons")[[1]]),1)),as.integer(1))
      proj4string(temp)<-CRS(proj4string(geoprj))
      shape<-as(temp,"owin")
    }
    xrng<-shape$xrange
    yrng<-shape$yrange
#    cat("\t",xrng,yrng,"\n")
#    plot(shape)
    hc<-1                                       #Count for Halton sequence
    hszv<-rep(pop[[i]][,1],times=pop[[i]][,2])  #Create household size vector
    nhh<-length(hszv)                           #Get number of households
    if(nhh>0){
      if(nhh>1)            
        hszv<-sample(hszv)                        #Scramble order
      hloc<-matrix(NA,nhh,3)                    #Household location matrix
      for(j in 1:nhh){                          #Walk households
        flag<-FALSE
        while(!flag){                           #Find a legal x,y loc
          temp<-switch(match.arg(method),
            uniform=runif(2,c(xrng[1],yrng[1]),c(xrng[2],yrng[2])),
            halton=rqhalton(1,d=2,min=c(xrng[1],yrng[1]), max=c(xrng[2],yrng[2]),base=c(2,3),seq.start=hc)
          )
          hc<-hc+1
          if(inside.owin(temp[1],temp[2],shape))
            flag<-TRUE
        }
        if(j>1){                                #Handle vertical stacking
          hd<-rowSums(sweep(hloc[1:(j-1),-3,drop=FALSE],2,temp,"-")^2)^0.5
          if(any(hd<stack.rad)){
            temp<-c(temp,max(0,hloc[1:(j-1),3][hd<stack.rad])+stack.dis)
            hloc[j,]<-temp
          }else{
            temp<-c(temp,0)
            hloc[j,]<-temp
          }
        }else{
          temp<-c(temp,0)
          hloc[j,]<-temp
        }
        for(k in 1:hszv[j]){   #Walk the household members
          flag<-FALSE
          while(!flag){                       #Find a legal x,y perturbation
            temp2<-runif(2,c(0,0),c(household.jitter,2*pi))
            temp2<-c(temp[1]+temp2[1]*sin(temp2[2]), temp[2]+temp2[1]*cos(temp2[2]))
            if(inside.owin(temp2[1],temp2[2],shape))
              flag<-TRUE
          }
          x[pc]<-temp2[1]
          y[pc]<-temp2[2]
          z[pc]<-temp[3]
          pc<-pc+1
        }
      } 
    }
  }
  #Return the results
  list(x=x,y=y,z=z)
}
