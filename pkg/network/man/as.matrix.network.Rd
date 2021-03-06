\name{as.matrix.network}
\alias{as.matrix.network}
\alias{as.matrix.network.adjacency}
\alias{as.matrix.network.incidence}
\alias{as.matrix.network.edgelist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Coerce a Network Object to Matrix Form }
\description{
  The \code{as.matrix} methods attempt to coerce their input to a matrix in adjacency, incidence, or edgelist form.  Edge values (from a stored attribute) may be used if present.
}
\usage{
\method{as.matrix}{network}(x, matrix.type = NULL, attrname = NULL, ...)
\method{as.matrix.network}{adjacency}(x, attrname=NULL, 
    expand.bipartite = FALSE, ...)
\method{as.matrix.network}{edgelist}(x, attrname=NULL, 
    as.sna.edgelist = FALSE, ...)
\method{as.matrix.network}{incidence}(x, attrname=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{network} }
  \item{matrix.type}{ one of \code{"adjacency"}, \code{"incidence"}, \code{"edgelist"}, or \code{NULL} }
  \item{attrname}{ optionally, the name of an edge attribute to use for edge values }
  \item{expand.bipartite}{ logical; if \code{x} is bipartite, should we return the full adjacency matrix (rather than the abbreviated, two-mode form)?}
  \item{as.sna.edgelist}{logical; should the edgelist be returned in sna edglist form?}
  \item{...}{ additional arguments. }
}
\details{
  If no matrix type is specified, \code{\link{which.matrix.type}} will be used to make an educated guess based on the shape of \code{x}.  Where edge values are not specified, a dichotomous matrix will be assumed.

  Edgelists returned by these methods are by default in a slightly different form from the \code{sna} edgelist standard, but do contain the \code{sna} extended matrix attributes (see \code{\link{as.network.matrix}}).  They should typically be compatible with \code{sna} library functions.  To ensure compatibility, the \code{as.sna.edgelist} argument can be set (which returns an exact \code{sna} edgelist). The \code{\link{as.edgelist}} function also returns a similar edgelist matrix but with an enforced sorting. 
  
  If the \code{attrname} attribute is used to include a charcter attribute, the resulting edgelist matrix will be character rather than numeric.
  
  Note that adjacency matrices may also be obtained using the extraction operator.  See the relevant man page for details.
}
\value{
An adjacency, incidence, or edgelist matrix
}
\references{ Butts, C. T.  (2008).  \dQuote{network: a Package for Managing Relational Data in R.}  \emph{Journal of Statistical Software}, 24(2).  \url{http://www.jstatsoft.org/v24/i02/}
}
\author{ Carter T. Butts \email{buttsc@uci.edu} and David Hunter \email{dhunter@stat.psu.edu} }
%\note{ ~~further notes~~ }
%
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{which.matrix.type}}, \code{\link{network}}, \code{\link{network.extraction}},\code{\link{as.edgelist}}}
\examples{
#Create a random network
m <- matrix(rbinom(25,1,0.5),5,5)
diag(m) <- 0
g <- network(m)

#Coerce to matrix form
as.matrix.network(g,matrix.type="adjacency")
as.matrix.network(g,matrix.type="incidence")
as.matrix.network(g,matrix.type="edgelist")

#Can also use the extraction operator
g[,]                            #Get entire adjacency matrix
g[1:5,6:10]                     #Obtain a submatrix
}
\keyword{ classes }% at least one, from doc/KEYWORDS
\keyword{ graphs }% __ONLY ONE__ keyword per line
