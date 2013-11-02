\name{NadaiPlot}
\alias{NadaiPlot}
\title{
Nadai plotting tool for fabric and strain data.
}
\description{
Function creates a standardized Nadai plot to better display calculated octahedral shear strain and Lode parameters of fabric and strain ellipsoids. 
}
\usage{
NadaiPlot(oss = 1, lp = 0, out.file,
  oss.int = 1, max.oss = 3,
  plot.title = "Nadai plot", labs)
}

\arguments{
  \item{oss}{
Octahedral shear strain.
}
  \item{lp}{
Lode parameter.
}
  \item{out.file}{
Path to PDF file for plotting.
}
  \item{oss.int}{
Interval used between arcs of equal octahedral shear strain in plot.
}
  \item{max.oss}{
Maximum octahedral shear strain to be plotted at top of graph (used to scale the plot). Value must be an integer.
}
  \item{plot.title}{
Title of the plot
}
  \item{labs}{
Vector of labels of length equal to the length of \code{oss} and \code{lp} to add to plot.
}
}

\value{
Creates PDF of Nadai plot.
}
\references{
Nadai, A. (1963) Theory of fracture and flow of solids.
}
\author{
Jeffrey R. Webber
}
