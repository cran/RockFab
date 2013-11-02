\name{ParEx}
\alias{ParEx}

\title{
Parameter extraction tool to obtain location, eccentricity, and orientation of objects from sectional image traces.
}
\description{
Function reads in an image file as accepted by the \code{readImage} function from the \code{EBImage} package and uses routines therein to obtain location, eccentricity, and orientation data.
}

\usage{
ParEx(img.file, out.file, section.title)
}

\arguments{
  \item{img.file}{
Path to image file of sectional object traces.
}
  \item{out.file}{
Path of PDF file to plot fitted object ellipses.
}
  \item{section.title}{
Title of PDF plot.
}
}

\value{
Function creates a PDF plot of fitted ellipses and returns a data frame object containing the parameters of each object.
}
\references{
Pau, G., A. Oles, M. Smith, O. Sklyar, and W. Huber. (2013) EBImage: Image processing toolbox for R.
}
\author{
Jeffrey R. Webber
}
\note{
All image analysis routines are obtained from the EBImage package
}
