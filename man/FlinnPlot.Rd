\name{FlinnPlot}
\alias{FlinnPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to plot a linear Flinn diagram from octahedral shear strain and Lode parameters of ellipsoids. 
}
\description{
Function produces a standardized linear Flinn diagram from the octahedral shear strains and Lode parameters of fitted ellipsoids. 
}
\usage{
FlinnPlot(oss = 1, lp = 0, out.file, max.k = 5,
plot.title = "Flinn diagram", labs)
}
\arguments{
  \item{oss}{
A vector of octahedral shear strain values.
}
  \item{lp}{
A vector of Lode parameters
}
  \item{out.file}{
Path to PDF file for plotting.
}
  \item{max.k}{
Maximum x and y plotting limits. Either X/Y or Y/Z. Default 5 is good for moderate strains. Higher strains may need larger values.
}
  \item{plot.title}{
Title of the plot.
}
  \item{labs}{Vector of labels of length equal to the length of \code{oss} and \code{lp} to add to plot.
}
}
\value{
Creates a PDF of Flinn diagram.
}
\references{
Flinn, Derek. "On the deformation of the Funzie conglomerate, Fetlar, Shetland." The Journal of Geology (1956): 480-505.
}
\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{NadaiPlot}}
}
