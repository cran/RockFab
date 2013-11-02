\name{FryData}
\alias{FryData}
\title{
Determine Fry points and scale to central void region.
}
\description{
Function takes locations of object centroids and calculates the parameters required to execute an interactive Fry plot. This function also determines a scale parameter that aids in data reduction.
}
\usage{
FryData(object.data, pie.step = 5, expansion = 2, section.name)
}

\arguments{
  \item{object.data}{
  Data frame object returned from \code{\link{ParEx}}.
}
  \item{pie.step}{
Angular distance in degrees as a multiple of 180. This is used to divide the fry plot into n wedge shaped regions of distance \code{pie.step} where the closest point to the origin within each wedge is saved and the largest of these values is returned. This provides an estimate of the central void diameter.
}
  \item{expansion}{
An expansion factor to be multiplied to the central void diameter estimate in order to properly scale Fry plots.
}
  \item{section.name}{
The name of the sectional analysis. 
}
}

\value{
Returns an object of class \code{\linkS4class{FRY}}.
}

\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{FryNormData}},\code{\link{RfPhi}}
}
