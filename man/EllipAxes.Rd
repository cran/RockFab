\name{EllipAxes}
\alias{EllipAxes}

\title{
Calculate ellipsoid axial lengths based on octahedral shear strain and Lode parameter.
}
\description{
Function uses the octahedral shear strain and Lode parameter of a desired strain ellipsoid and returns the normalized axial lengths X Y and Z.
}
\usage{
EllipAxes(es, nu)
}

\arguments{
  \item{es}{
Octahedral shear strain. Values must be positive.
}
  \item{nu}{
Lode parameter. Values must be between -1 and 1.
}
}
\value{
A numeric vector of length three with values returned in descending order (i.e. X, Y, and Z)
}
\references{
See for example: Ramsay, J. and M. Huber (1993). The techniques of modern structural geology.
}
\author{
Jeffrey R. Webber
}
\note{
Not used in RockFab scripts but can be useful for other endeavors.
}


\examples{
es <- runif(min = 0, max = 3, n = 1)
nu <- runif(min = -1, max = 1, n = 1)

EllipAxes(es = es, nu = nu)

}