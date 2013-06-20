\name{EDdata}
\alias{EDdata}
\docType{data}
\title{
EDdata is used for exemplifying the usage of function RadialPlotter
}
\description{
35 equivalent dose values come from Tengger Desert in northern china
}
\usage{data(EDdata)}
\format{
A data.frame with two columns contained equivalent dose and associated errors
}
\source{
OSL dated by Peng Jun in Chinese Academy of Sciences, Qinghai Insititute of Salt Lakes
}
\references{
Unpublished OSL equivalent doses data
}
\examples{
#
# loading the equivalent dose data
  data(EDdata)
#
# plot a probability density curve 
  x=EDdata
  plot(density(x[,1], weights=x[,2]/sum(x[,2])),col='blue',
       main='Sample.Tengger',xlab='De(Gy)',ylab='Density')
}
\keyword{package RadialPlotter}
\keyword{OSL dating}
\keyword{OSL age models}
