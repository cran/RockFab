\name{StereoPoint}
\alias{StereoPoint}
\title{
Add line data to stereonet.
}
\description{
Adds point to stereonet based on the trend and plunge of linear data.
}
\usage{
StereoPoint(my.az = 90, my.inc = 45,
  my.color = "black", my.pch = 19,
  my.size = 0.25, my.label)
}

\arguments{
  \item{my.az}{
Lineation trend as an azimuth.
}
  \item{my.inc}{
Lineation plunge or inclination.
}
  \item{my.color}{
Point color.
}
  \item{my.pch}{
pch value. See ?par
}
  \item{my.size}{
Size of points. Calls cex(my.size). See ?par
}
  \item{my.label}{
Label to add next to point
}
}
\value{
Adds point to stereonet based on the trend and plunge of linear data.
}
\references{
Modified from the package RFOC by Jonathan M. Lees.
}
\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{StereoPlot}}, \code{\link{StereoWeb}}, \code{\link{StereoPlane}}, \code{\link{StereoCirc}}
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
