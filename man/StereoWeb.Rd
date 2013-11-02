\name{StereoWeb}
\alias{StereoWeb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Equal area stereonet grid plotter.
}
\description{
Creates an equal-area net to existing stereographic projection.
}
\usage{
StereoWeb()
}

\value{
Plots an equal-area net to existing stereographic projection with 10 degree cells.
}
\references{
Modified from the package RFOC by Jonathan M. Lees.
}
\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{StereoPlot}}, \code{\link{StereoCirc}}, \code{\link{StereoPlane}}, \code{\link{StereoPoint}}
}
\examples{
#Plot random plane and lineation
strike <- runif(min = 0, max = 360, n = 1)
dip <- runif(min = 0, max = 90, n = 1)
trend <- runif(min = 0, max = 360, n = 1)
plunge <- runif(min = 0, max = 90, n = 1)

if(interactive()){
  StereoPlot(my.title = "Example Stereonet")
  StereoWeb()
  StereoCirc()
  StereoPoint(my.az = trend, my.inc = plunge)
  StereoPlane(my.strike = strike, my.dip = dip)
}
}
