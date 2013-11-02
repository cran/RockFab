\name{RFPHI-class}
\Rdversion{1.1}
\docType{class}
\alias{RFPHI-class}

\title{Class \code{"RFPHI"}}
\description{
Class to hold basic RfPhi data.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RFPHI", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{sectionName}:}{Object of class \code{"character"} defining the section ID string.}
    \item{\code{vectorMean}:}{Object of class \code{"numeric"} containing the calculated vector mean of object long axes.}
    \item{\code{harmonicMean}:}{Object of class \code{"numeric"} containing the harmonic mean of object axial ratios.}
    \item{\code{strainRatio}:}{Object of class \code{"numeric"} containing the calculated sectional strain or fabric ratio.}
    \item{\code{indexSymmetry}:}{Object of class \code{"numeric"} containing the index of symmetry based on the vector and harmonic means.}
    \item{\code{sampleSize}:}{Object of class \code{"numeric"} containing the number of objects used in the analysis. }
    \item{\code{meanObjectArea}:}{Object of class \code{"numeric"} containing the average object area based on fitted object ellipses.}
    \item{\code{rsAxes}:}{Object of class \code{"numeric"} of length two containing the axial lengths of the calculated sectional strain or fabric ellipse.}
    \item{\code{chiSquare}:}{Object of class \code{"data.frame"} containing the results of the chi squared statistical test for uniform distribution at each increment of anti-strain.}
  }
}
\section{Methods}{
No methods defined with class "RFPHI" in the signature.
}

\author{
Jeffrey R. Webber
}

\seealso{
\code{\linkS4class{FRY}}
}
\examples{
showClass("RFPHI")
}
\keyword{classes}
