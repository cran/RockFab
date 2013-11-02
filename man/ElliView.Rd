\name{ElliView}
\alias{ElliView}
\title{
Interactive 3D fabric and strain ellipsoid visualization tool.
}
\description{
Function deforms a unit sphere based on the parameters of the calculated fabric or strain ellipsoid and opens an rgl device to interactively view the result.
}
\usage{
ElliView(elli.dat, out.file, obj.res = 3)
}

\arguments{
  \item{elli.dat}{
Data frame object returned from \code{\link{ElliIn}}.
}
  \item{out.file}{
PLY file to save ellipsoid object to.
}
  \item{obj.res}{
Number of divisions to apply to spherical mesh. Larger numbers better approximate the ellipsoid surface but at the expense of file size and render time. 
}
}
\value{
PLY 3D object file.
}

\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{ElliIn}}
}
