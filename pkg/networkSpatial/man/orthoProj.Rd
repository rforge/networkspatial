\name{orthoProj}
\alias{orthoProj}
\title{Converts SpatialPolygons in Longituted/Latitude to an Orthogonal projection...}
\usage{orthoProj(polygon, pp)}
\description{Converts SpatialPolygons in Longituted/Latitude to an Orthogonal projection}
\details{Primarly for use with the \code{\link{gplot.spatial}} and \code{\link{rspop.SpatialPolygonDataFrame}}}
\value{\code{\link{SpatialPolygons}}}
\note{Should be used with UScensus* packages}
\references{See Almquist, Zack W. (2010). US Census Spatial and Demographic Data in R: The UScensus2000 Suite of Packages. Journal of Statistical Software, 37(6), 1--31.}
\author{Zack W. Almquist \email{TBD@umn.edu}}
\arguments{\item{polygon}{\code{\link{SpatialPolygons}}}
\item{pp}{\code{\link{vector}} Longituted/Latitude point for the orthogonal projection}}
\examples{\dontrun{

poly<-minnesota.cdp[minnesota.cdp$name=="Minneapolis",]
buildBB(poly)
}}
