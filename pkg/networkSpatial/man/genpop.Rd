\name{genpop}
\alias{genpop}
\title{Builds a list of household population for a polygon...}
\usage{genpop(spa, hh, total.pop, hh.size, verbose=FALSE)}
\description{Builds a list of household population for a polygon}
\details{Primarily for use with \code{\link{rnspatial}}}
\value{\code{list}}
\note{Should be used with UScensus* packages}
\references{See Almquist, Zack W. (2010). US Census Spatial and Demographic Data in R: The UScensus2000 Suite of Packages. Journal of Statistical Software, 37(6), 1--31.}
\author{Carter T. Butts \email{buttsc@uci.edu}, Zack W. Almquist \email{TBD@umn.edu}}
\arguments{\item{spa}{\code{\link{SpatialPolygons}}}
\item{hh}{\code{vector} of household variables in \code{\link{SpatialPolygons}}}
\item{total.pop}{\code{string} of total population variable in \code{\link{SpatialPolygons}}}
\item{hh.size}{\code{numeric} number of possible household size catagories}
\item{verbose}{default \code{FALSE}}}
\examples{\dontrun{require(UScensus2000tract);data(minnesota.tract);gp<-genpop(minnesota.tract,hh=c("hh.1person","hh.2person","hh.3person", "hh.4person","hh.5person","hh.6person"),total.pop="pop2000",hh.size=6)}}
