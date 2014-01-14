\name{FryObjective}
\alias{FryObjective}

\title{
Objective ellipse fitting to central void.
}
\description{
Function selects nearest points to Fry plot origin in angular wedges and fits an ellipse to these points using a deterministic least
squares method. The wedges are iteratively modified from the resultant ellipse to better define selection bins that are narrower near the apogee of the central void and wider near the perigee.
}
\usage{
FryObjective(object.data, n.pass = 15, pie.step = 12,
expansion = 1.5, pie.pts = 1, section.name, ave.piepts = FALSE,
norm = TRUE)
}

\arguments{
  \item{object.data}{
  Data frame object returned from \code{\link{ParEx}}.
}
  \item{n.pass}{
Number of ellipse fitting iterations.
}
  \item{pie.step}{
Angular distance in degrees as a multiple of 360. This is initially used to divide the fry plot into n wedge shaped regions of distance \code{pie.step} where the closest n \code{pie.pts} to the origin within each wedge is saved. Used in initial ellipse fit.
}
  \item{expansion}{
An expansion factor to be multiplied to the central void diameter estimate in order to properly scale Fry plots.
}
  \item{pie.pts}{
The number of nearest points to the plot origin that are selected in each wedge bin.
}
  \item{section.name}{
The name of the sectional analysis.
}
  \item{ave.piepts}{
Logical, if TRUE, the coordinates of n \code{pie.pts} selected from each wedge bin are averaged prior to least squares ellipse fitting.
}
  \item{norm}{
Logical, if TRUE, a normalized Fry method is used.
}
}
\details{
This function creates either a normalized or standard Fry plot data frame from the centroids of object traces and attempts to objectively fit an ellipse to the central void. This is done initially by determining the nearest n \code{pie.pts} to the origin within each angular bin of length \code{pie.step} in degrees. A deterministic least squares method of ellipse fitting (modified code written by Michael Bedward) is applied to these points in order to make an initial estimate of the central void. This initial estimate will likely underestimate the Rs ratio because a bias of more points near the perigees of the central void and fewer points near the apogees. The bins are modified based on the estimate ellipse so that so that a smaller angular distance is used near the apogees and a larger angular distance near the perigees. Points are selected from these bins are previously described and a new ellipse is fitted. This procedure is repeated \code{n.pass} times to return the parameters of the objectively fitted central void ellipse.
}
\value{
Returns an object of class \code{\linkS4class{FRY}}.
}
\references{
Deterministic least squares method of ellipse fitting function is modified from code written by Michael Bedward.

Halir, Radim, and Jan Flusser. "Numerically stable direct least squares fitting of ellipses." Proc. 6th International Conference in Central Europe on Computer Graphics and Visualization. WSCG. Vol. 98. 1998.

Fry, N. (1979) Random point distributions and strain measurement in rocks.

Erslev, E. (1988) Normalized center-to-center strain analysis of packed aggregates.
}
\author{
Jeffrey R. Webber
}

\seealso{
\code{\link{RfPhi}}, \code{\link{FryData}}, 
\code{\link{FryNormData}},\code{\link{FrySelect}}
}
