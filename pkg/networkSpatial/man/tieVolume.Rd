\name{tieVolume}
\alias{tieVolume}
\alias{tieVolume.SpatialPolygons}
\alias{tieVolume.SpatialPolygonsDataFrame}
\alias{tieVolume.default} 
\alias{tieVolumecore}
\alias{tieVolumertree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
tieVolume computation
}
\description{
Computes tieVolume
}
\usage{

\method{tieVolume}{default}(poly,...) 

\method{tieVolume}{SpatialPolygons}(poly, pop, poly2 = NULL, pop2 = NULL, param, model = c("powlaw", "atpowlaw", "explaw", "loglaw", "atanlaw", "tflaw", "fflaw", "cplaw"), tol.pa = c(0.01, 1e-08), tol.mc = c(0.01, 1e-08), trunc = 1e-10, gd.mode = c("ellipsoid", "sphere"), nq.mode = c("halton", "uniform"), maxiter = 1e+05, log = FALSE,method="core",...) 

\method{tieVolume}{SpatialPolygonsDataFrame}(poly, pop, poly2 = NULL, pop2 = NULL, param, model = c("powlaw", "atpowlaw", "explaw", "loglaw", "atanlaw", "tflaw", "fflaw", "cplaw"), tol.pa = c(0.01, 1e-08), tol.mc = c(0.01, 1e-08), trunc = 1e-10, gd.mode = c("ellipsoid", "sphere"), nq.mode = c("halton", "uniform"), maxiter = 1e+05, log = FALSE,method="core",...) 

tieVolumecore(poly, pop, poly2 = NULL, pop2 = NULL, param, model = c("powlaw", "atpowlaw", "explaw", "loglaw", "atanlaw", "tflaw", "fflaw", "cplaw"), tol.pa = c(0.01, 1e-08), tol.mc = c(0.01, 1e-08), trunc = 1e-10, gd.mode = c("ellipsoid", "sphere"), nq.mode = c("halton", "uniform"), maxiter = 1e+05, log = FALSE)

tieVolumertree(poly, pop, param, model = c("powlaw", "atpowlaw", "explaw", "loglaw", "atanlaw", "tflaw", "fflaw", "cplaw"), tol.pa = c(0.01, 1e-08), tol.mc = c(0.01,1e-08), trunc = 1e-10, gd.mode = c("ellipsoid", "sphere"), nq.mode = c("halton", "uniform"), maxiter =1e+05, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{poly}{\code{\link{SpatialPolygons}}}
\item{pop}{\code{list}}
\item{poly2 }{\code{\link{SpatialPolygons}}}
\item{pop2}{\code{list}}
\item{param}{\code{vector}}
\item{model}{\code{string}}
\item{tol.pa}{\code{numeric}}
\item{tol.mc}{\code{numeric}}
\item{trunc}{\code{numeric}}
\item{gd.mode}{\code{string}}
\item{nq.mode}{\code{string}}
\item{maxiter}{\code{numeric}}
\item{log}{\code{logical}}
\item{method}{\code{string}}
\item{...}{Other things}
}
\details{
\bold{Definitions}

Definition (Internal Tie Volume). Given a graph G with associated socio-physical coordinate space S and a region \eqn{A \subseteq S}, the internal tie volume of A with respect to G is defined as \eqn{V^i(A)=|E(G[A])|}.

Under the spatial Bernoulli graph model of Equation 2.3, it is possible to express the expectation of \eqn{V^i} thusly:

\eqn{E(V^i(A)) = E(|E(G[A])|) = \frac{1}{2} \int_A \int_A f_p(v_1)f_p(v_2)F_d(d(v_1,v_2))dv_1dv_2}

(where dv refers, as usual, to the differential operator; for consistency with Chapter 2, we retain the use of \eqn{d(v_i,v_j)} to refer to a distance function, on the assumption that context will render the meaning of the symbols unambiguous).


Definition (External Tie Volume). Given a graph G with associated socio-physical coordinate space S and a region \eqn{A \subseteq S}, the external tie volume of A with respect to G is defined as \eqn{V^e(A)=|E(G[A,\bar{A}])|}.

As with the internal tie volume, the expectation of this measure may be straightforwardly expressed as a sum of the expectations of tie variables in the usual manner:

\eqn{E(V^e(A)) = E(|E(G[A,\bar{A}])|) = \frac{1}{2} \int_A \int_\bar{A} f_p(v_1)f_p(v_2)F_d(d(v_1,v_2))dv_1dv_2}

Note that for population density functions which are sufficiently poorly behaved, the above may not be finite for infinite spaces. While this is clearly not a problem for empirical research (the world being a decidedly finite place), an exceptionally carelessly designed simulation could potentially be distorted by "edge effects" if populations tend to grow rapidly near the bounds of the simulated geography.

\bold{Spatial Bernoulli Graphs: param and model}

For details on the possible models \code{tieVolume} can handle see \code{\link{rnspatial}}.

\bold{rtree method}

}
\value{
Matrix object.
}
\references{
 Butts, Carter T. (forthcoming). Space and Structure: Methods and Models for Large-
Scale Interpersonal Networks. Springer, under contract.
}
\author{
Carter T. Butts, University of California, Irvine <buttsc@uci.edu>
}
\seealso{
See Also as \code{\link{rnspatial}}.
}
\examples{
\dontrun{
### Requires library(UScensus2000)
library(UScensus2000)
data(oregon.cdp)
bend<-city("bend","or")
bend.tv<-tieVolume(bend,bend$pop2000, param=c(0.937,0.538/1000,2.956),model="atpowlaw",tol.mc=c(.01,1e-2),tol.pa=c(.01,1e-2),nq.mode="uniform", gd.mode="sphere")

bend.tv.2<-tieVolume(bend,bend$pop2000, param=c(0.937,0.538/1000,2.956),poly2=bend[1,],pop2=1,model="atpowlaw",tol.mc=c(.01,1e-2),tol.pa=c(.01,1e-2),nq.mode="uniform", gd.mode="sphere")

bend.tv.r<-tieVolume(bend,bend$pop2000, param=c(0.937,0.538/1000,2.956),model="atpowlaw",tol.mc=c(.01,1e-2),tol.pa=c(.01,1e-2),nq.mode="uniform", gd.mode="sphere",method="rtree")

bend2<-as(bend,"SpatialPolygons")
bend.tv.sp<-tieVolume(bend2,rep(1,length(bend2@polygons)), param=c(0.937,0.538/1000,2.956),model="atpowlaw",tol.mc=c(.01,1e-2),tol.pa=c(.01,1e-2),nq.mode="uniform", gd.mode="sphere")
bend.tv.sp2<-tieVolume(bend2,rep(1,length(bend2@polygons)), param=c(0.937,0.538/1000,2.956),model="atpowlaw",tol.mc=c(.01,1e-2),tol.pa=c(.01,1e-2),nq.mode="uniform", gd.mode="sphere",method="rtree")

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
