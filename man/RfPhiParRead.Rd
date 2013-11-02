\name{RfPhiParRead}
\alias{RfPhiParRead}

\title{
Read manual data collected from the RfPhi method.
}
\description{
This function creates a data frame object typically constructed from sectional image data when manual Rf/Phi data is used. 
}
\usage{
RfPhiParRead(rfphi.file)
}

\arguments{
  \item{rfphi.file}{
A tab delineated text file containing the length of the long axis, short axis, and the orientation (rake) of objects measured in the field using a convention of rake measured positive clockwise from the direction of strike following right hand rule. Text file header must be used and should be defined as follows: l.axis  s.axis  phi
}
}

\value{
Returns a data frame object for use as input in \code{\link{RfPhi}}.
}

\author{
Jeffrey R. Webber
}


\seealso{
\code{\link{RfPhi}}
}
