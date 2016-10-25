\name{genpopIndv}
\alias{genpopIndv}
\title{Builds a list of indivuals for a polygon...}
\usage{genpopIndv(spa, total.pop)}
\description{Builds a list of indivuals for a polygon}
\details{Primarly for use with \code{\link{rnspatial}}}
\value{\code{list}}
\note{Should be used with UScensus* packages}
\references{See Almquist, Zack W. (2010). US Census Spatial and Demographic Data in R: The UScensus2000 Suite of Packages. Journal of Statistical Software, 37(6), 1--31.}
\author{Zack W. Almquist \email{TBD@umn.edu}}
\arguments{\item{spa}{\code{\link{SpatialPolygons}}}
\item{total.pop}{\code{string} of total population variable in \code{\link{SpatialPolygons}}}}
\examples{\dontrun{require(UScensus2000tract);data(minnesota.tract);gp<-genpopIndv(minnesota.tract,"pop2000")}}
