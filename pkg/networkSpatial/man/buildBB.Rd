\name{buildBB}
\alias{buildBB}
\title{Builds a bounding block polygon...}
\usage{buildBB(polygon, bb.epsilon=0)}
\description{Builds a bounding block polygon}
\details{Primarly for use with the \code{\link{gplot.spatial}}}
\value{\code{\link{SpatialPolygons}}}
\note{Should be used with UScensus* packages}
\references{See Almquist, Zack W. (2010). US Census Spatial and Demographic Data in R: The UScensus2000 Suite of Packages. Journal of Statistical Software, 37(6), 1--31.}
\author{Zack W. Almquist \email{TBD@umn.edu}}
\arguments{\item{polygon}{\code{\link{SpatialPolygons}}}
\item{bb.epsilon}{\code{numeric} default set to 0}}
\examples{\dontrun{

poly<-minnesota.cdp[minnesota.cdp$name=="Minneapolis",]
buildBB(poly)
}}
