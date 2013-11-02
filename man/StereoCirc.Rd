\name{StereoCirc}
\alias{StereoCirc}
\title{
Primitive circle plotter for stereonet.
}
\description{
This function simply adds a primitive circle to an existing plot and is used in coordination with other stereonet plotting tools. 
}
\usage{
StereoCirc(n.seg = 360)
}
\arguments{
  \item{n.seg}{
The number of segments to use in circle plot. Defaults to 360.
}
}
\value{
Adds primitive circle to open plot.
}
\references{
Modified from the package RFOC by Jonathan M. Lees.
}
\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{StereoPlot}}, \code{\link{StereoWeb}}, \code{\link{StereoPlane}}, \code{\link{StereoPoint}}
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
