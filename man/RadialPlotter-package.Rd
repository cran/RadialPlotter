\name{RadialPlotter-package}
\alias{RadialPlotter-package}
\docType{package}
\title{Numeric routines for optical stimulated luminescence dating}
\description{
Five functions that wrapped with Fortran subroutines are available in package "RadialPlotter" currently.\cr\cr
Function calED() fits dose-response curve and calculates equivalent dose in optically stimulated luminescence dating, assiociated standard error of equivlalent dose can also be estimated using either simple transformation or monte Carlo method (Duller, 2007).\cr\cr
Function decomp() decomposes the OSL decay curve (type "CW" or "LM") to a specified number of components, which supports the combination of differential evolution and
Levenberg-Marquardt algorithm suggested by Bluszcz and Adamiec (2006).\cr\cr
Function fastED() approximates a fast-component equivalent dose with data abtained by the single aliquot regenerative (SAR) protocol (Murray and Wintle, 2000).\cr\cr
Function RadialPlotter() performs statistical age models (log-scale) analysis (Galbraith and Roberts, 2012) and radial plot drawing (Galbraith, 1988). Routines for statistical age models analysis include: \cr
Central age model (CAM) \cr
Finite mixture age model (FMM) \cr
Minimum age model with three parameters (MAM3) \cr
Minimum age model with four parameters (MAM4) \cr\cr
Function sgcED() is a wrapped version of function calED(), fitting OSL dose-response curve and calculating a series of equivalent doses using standardised growth curves (SGC) method (Roberts and Duller, 2004).
}
\details{
\tabular{ll}{
Package: \tab RadialPlotter\cr
Type: \tab Package\cr
Version: \tab 2.2\cr
Date: \tab 2013-11-27\cr
License: \tab GPL-3\cr
}
}
\author{
\tabular{ll}{
Peng Jun \tab University of Chinese Academy of Sciences, Beijing, China
}
\bold{Package maintainer}

Peng Jun
\email{pengjun10@mails.ucas.ac.cn}
}

\keyword{package RadialPlotter}
\keyword{OSL dating}
\keyword{OSL age models}
\keyword{decay curve decomposition}
\keyword{SAR method}
\keyword{SGC method}
