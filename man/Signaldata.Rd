\name{Signaldata}
\alias{Signaldata}
\docType{data}
\title{
Signaldata is used for exemplifying the usage of function decomp
}
\description{
OSL signal data of type "CW" and "LM".
}
\usage{data(Signaldata)}
\format{
  A list contains both OSL decay curve data of type "CW" and "LM"
  \describe{
    \item{cw}{a data.frame of OSL signal decay data of type "CW" measured over 40 s in 250 equal width channels, which is
the nature OSL signal of a sample from Tengger Desert in northern china (Peng and Han, 2013)}
    \item{lm}{a data.frame of synthesized OSL signal decay data of type "LM" over 3600s with a1=1e3, a2=1e3, a3=1e4, a4=1e4, a5=1e5, b1=10.5, b2=2.5, b3=0.5, b4=0.05, b5=0.005}
}
}
\references{
Peng, J., Han, F.Q., 2013. Selections of fast-component OSL signal using sediments from the south edge of Tengger Desert. Acta Geoscientica Sinica, 34 (1), pp. 1-7.
}
\examples{
data(Signaldata)
names(Signaldata)
}
\keyword{package RadialPlotter}
\keyword{OSL dating}
\keyword{CW-OSL decay curve}
\keyword{LM-OSL decay curve}
