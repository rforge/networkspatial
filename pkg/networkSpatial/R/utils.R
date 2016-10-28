######################################################################
#
# utils.R
#
# Written by Carter T. Butts <buttsc@uci.edu>.
#
# Last Modified 07/31/09
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/networkSpatial package
#
# This file contains various utilities needed for the networkSpatial
# package.
#
# Contents:
#  bbAreaWGS84
#  geoDistWGS84
#  listPrimes
#
######################################################################


#Calculate bounding box areas on the geosphere (or rectilinear coordinates)
bbAreaWGS84<-function(coord,mode=c("sphere","euclidean")){
  area<-as.double(0)
  dim<-as.integer(NCOL(coord))
  mode<-as.integer(switch(match.arg(mode),
    sphere=0,
    euclidean=1
  ))
  .C("bbAreaCalc_R",area=area,as.double(coord),dim,mode, PACKAGE="networkSpatial")$area;
}


#Calculate geodesic distances on the WGS84 ellipsoid (or sphere, if desired)
geoDistWGS84<-function(lat0,lon0,lat1,lon1,mode=c("ellipsoid","sphere")){
  n<-length(lat0)
  mode<-switch(match.arg(mode),
    ellipsoid=0,
    sphere=1
  )
  .C("geoDistWGS84_R",as.integer(n),gd=as.double(rep(0,n)),as.double(lat0), as.double(lon0),as.double(lat1),as.double(lon1),as.integer(mode), PACKAGE="networkSpatial")$gd
}


#List the first n primes that are >=start
listPrimes<-function(n,start=1){
  #require(spatstat)
  #Start with a precomputed list of the first 100 primes
  firstprimes<-c(1,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71, 73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523)
  #Find the first n primes >=start
  firstprimes<-firstprimes[firstprimes>=start]
  if(length(firstprimes)>0){
    primes<-firstprimes[1:min(length(firstprimes),n)]
    n<-n-length(primes)
  }
  if(n>0){   #Must have run out of primes
    if(length(primes)>0)
      start<-primes[length(primes)]+1
    if(start%%2==0)                       #Avoid evens
     start<-start+1
    while(n>0){                           #Hunt 'em the hard way
      if(is.prime(start)){
        primes<-c(primes,start)
        n<-n-1
      }
      start<-start+2                      #Can advance by 2, to avoid evens
    }
  }
  #Return the result
  primes
}
