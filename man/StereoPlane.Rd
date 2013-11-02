\name{StereoPlane}
\alias{StereoPlane}
\title{
Add plane data to stereonet.
}
\description{
Function adds a great circle to existing plot based on the strike and dip of a plane.
}
\usage{
StereoPlane(my.strike = 0, my.dip = 45, my.color = "black")
}
\arguments{
  \item{my.strike}{
The strike of the plane to plot as a great circle.
}
  \item{my.dip}{
The dip of the plane to plot as a great circle.
}
  \item{my.color}{
The line color used in the plot.
}
}

\value{
Adds great circle of plane to existing plot.
}
\references{
Modified from the package RFOC by Jonathan M. Lees.
}
\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{StereoPlot}}, \code{\link{StereoWeb}}, \code{\link{StereoCirc}}, \code{\link{StereoPoint}}
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
