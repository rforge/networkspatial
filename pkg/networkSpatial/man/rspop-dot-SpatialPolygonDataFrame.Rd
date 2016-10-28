\name{rspop.SpatialPolygonDataFrame}
\alias{rspop.SpatialPolygonDataFrame}
\alias{rspop}
\alias{rqhalton}
\title{Builds a random set of households or indivuals...}
\usage{rspop.SpatialPolygonDataFrame(geo, method=c("uniform",
    "halton"), stack.rad=10, stack.dis=4, household.jitter=5,
    longLat=TRUE, household=TRUE, list.pop)}
\description{Builds a random set of households or indivuals}
\details{Primarly for use with the \code{\link{rnspatial}}}
\value{\code{list}}
\note{Should be used with UScensus* packages}
\references{
Almquist, Zack W. (2010). US Census Spatial and Demographic Data in R: The UScensus2000 Suite of Packages. Journal of Statistical Software, 37(6), 1--31.

Almquist, Zack W. and Butts, Carter T. (2012). Point process models for household distributions within small areal units. Demographic Research, 26(12), 593--632.
}
\author{Zack W. Almquist \email{almquist@umn.edu}, Carter T. Butts \email{buttsc@uci.edu}}
\arguments{\item{geo}{\code{\link{SpatialPolygons}}}
\item{method}{\code{string}, "uniform" or "halton"}
\item{stack.rad}{\code{numeric} default set to 10}
\item{stack.dis}{\code{numeric} default set to 4}
\item{household.jitter}{\code{numeric} default set to 5}
\item{longLat}{\code{logical} default set to \code{TRUE}}
\item{household}{\code{logical} default set to \code{TRUE}}
\item{list.pop}{\code{list} requires necessary inputs for \code{\link{genpop}} and \code{\link{genpopIndv}}}}
\examples{\dontrun{

poly<-minnesota.cdp[minnesota.cdp$name=="Minneapolis",]
"hh.4person","hh.5person","hh.6person"),total.pop="pop2000",hh.size=6))

}}
