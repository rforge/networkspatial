\name{gplot.spatial}
\alias{gplot.spatial}
\title{Plot spatial embedded networks (undirected)...}
\usage{
gplot.spatial(x,coord,bb,list.edges,list.vertex,background=NULL,list.background,list.axis=list(one=list(side=1,cex.axis=.4,lwd=.5,las=1),two=list(side=2,cex.axis=.4,lwd=.5)),axis=TRUE,bg=NULL)
}
\description{Plot spatial embedded networks (undirected)}
\details{Primarily for use \code{\link{rnspatial}}}
\value{\code{NULL}}
\note{Should be used with UScensus* packages}
\references{See Almquist, Zack W. (2010). US Census Spatial and Demographic Data in R: The UScensus2000 Suite of Packages. Journal of Statistical Software, 37(6), 1--31.}
\author{Zack W. Almquist \email{TBD@umn.edu}}

\arguments{
\item{x}{\code{\link{network}}}
\item{coord}{\code{matrix}}
\item{bb}{\code{\link{SpatialPolygons}}}
\item{list.edges}{\code{list}}
\item{list.vertex}{\code{list}}
\item{background}{\code{\link{SpatialPolygons}}}
\item{list.background}{\code{list}}
\item{list.axis}{hi}
\item{axis}{\code{logical}}
\item{bg}{color}
}

\examples{\dontrun{

 require(UScensus2000tract)
 require(rgdal)
data(minnesota.tract)
data(minnesota.cdp)
minneapolis.cdp<-minnesota.cdp[minnesota.cdp$name=="Minneapolis",]
 index<-overlay(SpatialPoints(coordinates(minnesota.tract)),minneapolis.cdp)
minneapolis.tract<-minnesota.tract[which(!is.na(index)),]
coord<-rspop.SpatialPolygonDataFrame(geo=minneapolis.tract,method="uniform",list.pop=list(hh=c("hh.1person","hh.2person","hh.3person", "hh.4person","hh.5person","hh.6person"), total.pop="pop2000",hh.size=6),longLat=FALSE)
net<-rnspatial(1,cbind(coord$x,coord$y),c(0.937, 0.538, 2.956),model="atpowlaw")
projection.point<-coordinates(buildBB(minneapolis.tract))
coord<-SpatialPoints(cbind(coord$x,coord$y))
proj4string(coord)<-
 CRS(paste("+proj=ortho +lat_0=",projection.point[2], " +lon_0=", projection.point[1],collapse = "", sep = ""))
mpls.prj<-spTransform(minneapolis.tract,CRS(proj4string(coord)))
deg<-degree(net)
dego<-order(deg,decreasing =TRUE)
points<-coord[which(deg[dego[1000]]==deg)[1],]
index<-overlay(points,mpls.prj)
index<-unique(index)
index<-index[!is.na(index)]
poly<-buildBB(minneapolis.tract[index,])
poly<-spTransform(poly,CRS(proj4string(coord)))
net.coord<-coordinates(coord)
col<-rgb(1,0,0,.05)
col2<-rgb(1,0,0,.1)
png("mpls.png",width=1600,height=1200,res=300,bg="white")
par(mar=c(1, 1, 0, 0) + 0.1,mgp=c(0, 0.2,0),tcl=-.2)
gplot.spatial(net,net.coord,bb=bbox(mpls.prj), list.vertex=list(col=col,pch=19,cex=.01),list.edges=list(col=rgb(0,0,0,.3),lwd=.3),background=mpls.prj,list.background=list(border=rgb(0,0,1,.1)))
par(fig=c(.25, 4, .6, 3.1)/10,new=TRUE,mgp=c(0,0.2,0),tcl=-.2)
gplot.spatial(net,net.coord,background=mpls.prj,bb=bbox(poly), list.vertex=list(col=col2,pch=19,cex=.1), list.edges=list(col=rgb(0,0,0,.3)), list.background=list(border=rgb(0,0,1,.1)),bg="white")
dev.off()
}}
