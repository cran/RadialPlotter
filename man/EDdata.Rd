\name{EDdata}
\alias{EDdata}
\docType{data}
\title{
EDdata is used for exemplifying the usage of RadialPlotter
}
\description{
35 equivalent dose values come from Tengger desert in northern china
}
\usage{data(EDdata)}
\format{
A data.frame with two columns contained equivalent dose and associated errors
}
\source{
OSL dating by Peng Jun in Chinese Academy of Sciences, Qinghai Insititue of Salt Lakes
}
\references{
Unpublished OSL equivalent doses data
}
\examples{
#
# loading the equivalent dose data
  data(EDdata)
#
# plot a histogram
  hist(EDdata[,1])
}
\keyword{datasets}
