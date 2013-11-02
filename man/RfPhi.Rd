\name{RfPhi}
\alias{RfPhi}

\title{
Determine Rs and corresponding parameters following the RfPhi method
}
\description{
This function applies the RfPhi technique to object data following a Chi squared statistical test for the most uniform distribution of object orientations. 
}
\usage{
RfPhi(my.par, out.file, section.title, weight.vec = TRUE, unit.area = "px")
}

\arguments{
  \item{my.par}{
A data frame object returned by \code{\link{ParEx}} or \code{\link{RfPhiParRead}}.
}
  \item{out.file}{
Path to PDF file where a standardized plot of the sectional results is made.
}
  \item{section.title}{
Title of the sectional analysis. 
}
  \item{weight.vec}{
Logical, if TRUE an eccentricity weighting is applied such that more elliptical objects influence the vector mean stronger than objects with more circular shapes. 
}
  \item{unit.area}{
The unit of measure in which the mean object area is calculated. If, for example, the input data is returned from \code{\link{ParEx}} the unit is in pixels; however, if manual data is used, the unit used to measure axial lengths should be used.
}
}

\value{
Function creates a standardized PDF file containing an RfPhi plot, sectional Rs ellipse, Chi square plot, and a list of calculated properties. These properties are also returned via an object of class \code{\linkS4class{RFPHI}}.
}
\references{
Dunnet, D. (1969). A technique of finite strain analysis using elliptical particles.

Ramsay, J. (1969). Folding and fracturing of rocks.

Lisle, R. (1985). Geological strain analysis, a manual for the rf/phi method.

Mulchrone, K. and P. Meere (2001). A windows program for the analysis of tectonic strain using deformed elliptical markers.

Mookerjee, M. and S. Nickleach (2011). Three-dimensional strain analysis using Mathematica.
}
\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{FryData}}, \code{\link{FryNormData}}, \code{\link{FrySelect}}
}
