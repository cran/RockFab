\name{FRY-class}
\Rdversion{1.1}
\docType{class}
\alias{FRY-class}

\title{Class \code{"FRY"}}
\description{
A class to hold data from sectional Fry analyses.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FRY", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{sectionName}:}{Object of class \code{"character"} defining the section ID string.}
    \item{\code{vectorMean}:}{Object of class \code{"numeric"} for the calculated calculated central void rake.}
    \item{\code{strainRatio}:}{Object of class \code{"numeric"} for calculated strain ratio.}
    \item{\code{sampleSize}:}{Object of class \code{"numeric"} containing the number of objects used in sectional analysis.}
    \item{\code{rsAxes}:}{Object of class \code{"numeric"} with length of two containing the major and minor axial lengths of central void.}
    \item{\code{meanObjectArea}:}{Object of class \code{"numeric"} the average object area.}
    \item{\code{fryParams}:}{Object of class \code{"data.frame"} containing fry point coordinates, angle from origin, and distance from origin.}
    \item{\code{voidScale}:}{Object of class \code{"numeric"} containing the distance from origin plot should be scaled.}
  }
}
\section{Methods}{
No methods defined with class "FRY" in the signature.
}

\author{
Jeffrey R. Webber
}

\seealso{

\code{\linkS4class{RFPHI}}
}
\examples{
showClass("FRY")
}
\keyword{classes}
