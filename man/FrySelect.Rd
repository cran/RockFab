\name{FrySelect}
\alias{FrySelect}

\title{
Interactive central void selection tool for the Fry and normalized Fry techniques.
}
\description{
Creates an interactive Fry or normalized plot that enables a user to select the parameters of the central void.
}
\usage{
FrySelect(fry.data, max.dim, out.file, normalized = FALSE, select = TRUE)
}

\arguments{
  \item{fry.data}{
  An object of class \code{\linkS4class{FRY}} returned from \code{\link{FryData}} or \code{\link{FryNormData}}.
}
  \item{max.dim}{
Fry plot diameter.
}
  \item{out.file}{
Path to desired PDF file for results.
}
  \item{normalized}{
Logical, if TRUE, a normalized Fry plot is made.
}
  \item{select}{
Logical, if TRUE, manual selection is required. If FALSE, an object of \code{\linkS4class{FRY}} will be plotted.
}
}

\value{
Function produces a standardized PDF file containing the results of the sectional analysis and returns an object of class \code{\linkS4class{FRY}} with all parameters of the sectional analyses.
}
\references{
Fry, N. (1979) Random point distributions and strain measurement in rocks.

Erslev, E. (1988) Normalized center-to-center strain analysis of packed aggregates.
}
\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{RfPhi}}
}
