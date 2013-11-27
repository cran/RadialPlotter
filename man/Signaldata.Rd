\name{Signaldata}
\alias{Signaldata}
\docType{data}
\title{
Signaldata is used for exemplifying the usage of function decomp() and fastED()
}
\description{
OSL signal data of type "CW" and "LM".
}
\usage{data(Signaldata)}
\format{
  A list contains both OSL decay curve data of type "CW" and "LM"
  \describe{
    \item{cw}{a series of OSL decay curve of type "CW" obatined with SAR method that from Tengger Desert in northern china (Peng and Han, 2013)}
    \item{lm}{a data.frame of OSL signal decay data of type "LM" from Li and Li (2006)}
}
}
\references{
Li, S.H., Li, B., 2006. Dose measurement using the fast component of LM-OSL signals from quartz. Radiation Measurements, 41, pp.534-541.

Peng, J., Han, F.Q., 2013. Selections of fast-component OSL signal using sediments from the south edge of Tengger Desert. Acta Geoscientica Sinica, 34 (6), pp. 757-762.
}
\examples{
data(Signaldata)
names(Signaldata)
}
\keyword{package RadialPlotter}
\keyword{OSL dating}
\keyword{CW-OSL decay curve}
\keyword{LM-OSL decay curve}
