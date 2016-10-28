\name{rnspatial}
\alias{rnspatial}
\alias{bbAreaWGS84}
\alias{geoDistWGS84}
\alias{listPrimes}
\alias{rnspatial.grid} 
\alias{rnspatial.rtree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate Draws from a Spatial Bernoulli Graph Process }
\description{
  Simulates one or more draws from a spatial Bernoulli graph process, using one of several models.
}
\usage{
rnspatial(n, coord, param, accel.method = c("grid", "rtree"), ... )

rnspatial.grid(n, coord, param, model = c("powlaw", "atpowlaw", "explaw", 
    "loglaw", "atanlaw", "tflaw", "fflaw", "cplaw"), directed = FALSE,
    return.as.network = TRUE, minkowski.exp = 2, height.dimension = NULL,
    height.tol = 0, covar = NULL, psi = NULL, tol = 0.99, 
    epdiff.tol = 1e-09, subdivisions = 10, ... )

rnspatial.rtree(n, coord, param, model = c("powlaw", "atpowlaw", "explaw", 
    "loglaw", "atanlaw", "tflaw", "fflaw", "cplaw"), directed = FALSE,
    return.as.network = TRUE, minkowski.exp = 2, height.dimension = NULL,
    height.tol = 0, covar = NULL, psi = NULL, tol = 0.99, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ the number of draws to take. }
  \item{coord}{ a vector or matrix containing vertex coordinates (one row per vertex). }
  \item{param}{ a vector of parameters for the Bernoulli model. }
  \item{accel.method}{ the acceleration method to use when simulating the graph structure (currently gridding or R-trees).}
  \item{model}{ the edge probability model to use. }
  \item{directed}{ logical; should the generated graph be directed? }
  \item{return.as.network}{ logical; should the generated graph be returned as a network object?  (If \code{FALSE}, an \code{sna} edgelist is used.)}
  \item{minkowski.exp}{ the exponent to use for the Minkowski metric on vertex coordinates. }
  \item{height.dimension}{ optionally, the index of a dimension to interpret as a \dQuote{structure height} or artificial elevation when calculating distance.}
  \item{height.tol}{ for artificial elevations, how far apart on the non-elevation dimension can two vertices be before requiring a ground-level transit?}
  \item{covar}{ an optional matrix of vertex covariates to use when determining base tie probabilities at distance 0.}
  \item{psi}{ a vector of (logit) parameters to use with \code{covar}.}
  \item{tol}{ only bother to calculate tie probabilities (and draw edges) between two cells when the probability that there are no inter-cell ties is less than \code{tol}; higher numbers are more conservative, while lower numbers can (greatly) reduce simulation time. }
  \item{epdiff.tol}{ maximum edge probability difference to tolerate for use of the cell centroid approximation versus exact distances during simulation. }
  \item{subdivisions}{ the number of cells into which the vertex coordinate space should be separated for computational purposes (on each dimension). }
  \item{...}{ additional arguments.}
}
\details{

\code{rnspatial} can draw spatial Bernoulli graphs with a range of spatial interaction functions (selected by \code{model}).  Currently supported variants are as follows:

\describe{
\item{\code{powlaw}}{Power law model of the form }%\eqn{\frac{\theta_1}{(1+\theta_2 d)^{\theta_3}}}{theta_1/(1+theta_2 d)^theta_3}.}
\item{\code{atpowlaw}}{Attenuated power law of the form}% \eqn{\frac{\theta_1}{1+(\theta_2 d)^{\theta_3}}}{theta_1/(1+(theta_2 d)^theta_3)}.}
\item{\code{explaw}}{Exponential decay law of the form}% \eqn{\theta_1\exp(-\theta_2 d)}{theta_1 exp(-theta_2 d)}.}
\item{\code{loglaw}}{Logistic decay law of the form }%\eqn{\frac{(1+\theta_3)\theta_1}{1+\theta_3\exp(\theta_2 d)}}{(1 + theta_3) theta_1/(1 + theta_3 exp(theta_2 d))}.}
\item{\code{atanlaw}}{Arctangent law of the form }%\eqn{\theta_1(1-2/\pi \atan(theta_2d))}{theta_1 (1 - 2/Pi atan(theta_2 d))}.}
\item{\code{tfaw}}{Truncated linear law of the form}% \eqn{\max(\min(\theta_3 - \theta_2 d, \theta_1), \theta_4\theta1)}.}
\item{\code{fflaw}}{Fixed, \dQuote{disk-like} law in which tie probability is equal to}% \eqn{\theta_1} for \eqn{d<\theta_2}, otherwise \eqn{\theta_3 \theta1}.}
\item{\code{cplaw}}{Constant probability law (null model), in which tie probability is fixed at}% \eqn{\theta_1} for all \eqn{d}.}
}
\eqn{\theta}{theta} values in the above correspond to the entries of \code{param}, with \eqn{d} corresponding to the distance between vertices.  (In the directed case, both edges of the dyad are drawn independently.)  The distance employed can be any Minkowski metric on \code{coord}, the exponent being controlled by \code{minkowski.exp}.  

In addition to spatial effects, it is possible to scale the base tie probability between vertices by covariate effects.  If \code{covar} and \code{psi} are given, the spatially determined tie probability described above is multiplied by the inverse logit of \eqn{\psi |x_i-x_j|}, where \eqn{x_i,x_j} are vectors of covariates.  This provides a simple (albeit crude) means of adding non-spatial homophily effects and the like. 

\code{rnspatial} uses various tricks to accelerate the process of drawing networks on large coordinate sets.  The general method used is controlled by \code{accel.method}, which determines the specialized routine to be used when simulating the graph.  (One may also call these methods directly, via the appropriate variant function.)  For the \code{grid} method, the coordinate space is subdivided into a series of \dQuote{cells,} which are used to avoid explicit probability calculations for dyads which are sufficiently unlikely to contain edges.  The \code{subdivisions} argument controls the number of cells, with the space being split by into the indicated number of segments on each dimension.  (The total number of cells is thus \code{subdivisions^d}, where \code{d} is the dimensionality of \code{coord}.)  Computation is only performed for dyads spanning two cells when the probability that at least one cross-cell edge exists is no greater than \code{1-tol} -- adjusting \code{tol} thereby allows the user to balance the requirements of speed (low \code{tol}) and accuracy (high \code{tol}).  Where the difference in baseline probabilities across grid cells is less than \code{epdiff.tol}, then cell centroids are used instead of individual distances to determine tie probability (greatly reducing execution time, at cost of accuracy).  The ideal number of subdivisions to use varies with edge model and population distribution, and experimentation may be required to find a value which works well for a particular problem.  High values of \code{subdivision} will increase the precision with which avoidable dyads are identified, but will correspondingly increase the number of cell pairs which must be examined.  Use of the diagnostic timing information displayed during function execution is recommended. 

For \code{accel.method=="rtree"}, an R-tree traversal method is employed (using multidimensional bounding boxes).  The parameters for this method are as for \code{grid} above, except that the division of the coordinate set is fixed by the construction algorithm (and hence it has no analog to the \code{subdivisions} argument).  The spatial R-tree carries some initial overhead (construction is generally \eqn{O(N \log N)}{O(N log N)}, unlike the \eqn{O(N)} construction for the coordinate grid), but is adaptive; it may outperform the gridding method for spatially heterogeneous layouts, in high-dimensional spaces, etc.  As with the gridding method, a simple diagnostic of \dQuote{avoided} dyads is printed (although direct examination of timing is a better performance indicator).
}
\value{
  A network object containing the simulated graph, or a list of such objects if \code{n>1}.  If \code{return.as.network} is \code{FALSE}, \code{sna} edgelists are used instead.
}
\references{ Butts, Carter T.  (forthcoming).  Space and Structure: Models and Methods for Large-Scale Interpersonal Networks.  Springer.   }
\author{ Carter T. Butts \email{buttsc@uci.edu} }
\note{ Passing illegal parameter values in \code{param} could cause exciting behavior.  Caveat emptor.
}
\seealso{ \code{\link[network]{network}} }
\examples{
#Draw some 2D coordinates
coord<-matrix(runif(400,0,2),nc=2)

#Produce a graph using the attenuated power law
net<-rnspatial(1,coord,c(1,8,3),model="atpowlaw")

#Show the results
plot(net,coord=coord,vertex.cex=0.75)
}
\keyword{ graphs }% at least one, from doc/KEYWORDS
\keyword{ datagen }% __ONLY ONE__ keyword per line
