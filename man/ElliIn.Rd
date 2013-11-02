\name{ElliIn}
\alias{ElliIn}
\title{Read output files generated from the program Ellipsoid 2003.
}
\description{
ElliIn reads output files generated from the program Ellipsoid 2003 in order to create stereonets, Nadai plots, and interactive 3D ellipsoid models.
}
\usage{
ElliIn(elli.files)
}

\arguments{
  \item{elli.files}{Path of Ellipsoid 2003 file (.elli) to be read into R.
}
}
\value{
Returns a data frame object containing properties of fitted ellipsoids including the filename of the Ellipsoid 2003 file; X, Y, and Z axial lengths of the ellipsoids; X, Y, Z axial azimuths and inclinations; strike and dip of the XY plane; rake of lineation within XY plane; number of sections used to fit ellipsoid; incompatibility index; octahedral shear strain, and Lode parameter.
}
\references{
Launeau, P. and P. Robin (2005). Determination of fabric and strain ellipsoids from measured sectional ellipses implementation and applications. Journal of structural geology 27 (12), 2223--2233
}
\author{
Jeffrey R. Webber
}
\seealso{
\code{\link{ElliOut}}
}
