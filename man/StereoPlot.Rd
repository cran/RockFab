\name{StereoPlot}
\alias{StereoPlot}
\title{
Setup plot extent for stereonet.
}
\description{
Creates a standardized plot with predetermined extents and margins to used with other stereonet tools. 
}
\usage{
StereoPlot(my.title = "Stereonet", new = TRUE, pdf.file)
}

\arguments{
  \item{my.title}{
Title of the stereonet plot.
}
  \item{new}{
Logical, if FALSE, new plot will open in active window.
}
  \item{pdf.file}{
The filename to plot stereonet as a pdf. If missing, windows() is used. 
}
}

\value{
Sets up a plot for stereographic projections.
}
\references{
Modified from the package RFOC by Jonathan M. Lees.
}
\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{StereoCirc}}, \code{\link{StereoWeb}}, \code{\link{StereoPlane}}, \code{\link{StereoPoint}}
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
