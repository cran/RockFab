\name{ElliOut}
\alias{ElliOut}
\title{
Write results of sectional analyses to Ellipsoid 2003 compatible files.
}
\description{
ElliOut writes an ASCII file formated for use in the program Ellipsoid 2003 for ellipsoid fitting to sectional data.
}
\usage{
ElliOut(my.results, setup.file, out.file, use.ratio = FALSE)
}

\arguments{
  \item{my.results}{
An object or list of objects of class RFPHI or FRY.
}
  \item{setup.file}{
Location of setup.txt file containing orientation data corresponding to the sectional fabric or strain data. 
}
  \item{out.file}{
File to write Ellipsoid 2003 file to.
}
  \item{use.ratio}{
Logical, if TRUE, Rs ratio is used instead of axial lengths. Best used when axial lengths between sectional analyses are not meaningful. 
}
}

\value{
File for use in Ellipsoid 2003. 
}
\references{
Launeau, P. and P. Robin (2005). Determination of fabric and strain ellipsoids from measured sectional ellipses implementation and applications. Journal of structural geology 27 (12), 2223--2233
}
\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{ElliIn}}
}
