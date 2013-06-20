\name{CW.Signal.data}
\alias{CW.Signal.data}
\docType{data}
\title{
CW.Signal.data is used for exemplifying the usage of function decomp
}
\description{
CW OSL shine down signal data measured over 40 s in 250 equal width channels
}
\usage{data(CW.Signal.data)}
\format{
  A data frame with 250 observations on the following 2 variables.
  \describe{
    \item{\code{time.s.}}{a numeric vector}
    \item{\code{signal.count.}}{a numeric vector}
  }
}
\source{
OSL Signal data of a sample from Tengger Desert in northern china
}
\references{
Unpublished OSL data
}
\examples{
data(CW.Signal.data)
plot(CW.Signal.data,type='o')
}
\keyword{package RadialPlotter}
\keyword{OSL dating}
\keyword{CW decay curve}
