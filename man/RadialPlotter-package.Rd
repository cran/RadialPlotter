\name{RadialPlotter-package}
\alias{RadialPlotter-package}
\docType{package}
\title{Package for OSL equivalent dose analysis and radial plot drawing }
\description{
Three functions are availabe in package 'RadialPlotter' currently.\cr\cr
Function "decomp" decomposes the OSL decay curve to a specified number of components using a combination of differential evolution and Levenberg-Marquardt algorithm.\cr\cr
Function "RadialPlotter" performs statistical age models (log-scale) analysis and radial plot drawing. Routines for statistical age models analysis using fortran subroutines wrapped function 'RadialPlotter'. Facilitating the analysis of Galbraith's statistical age models, that are: \cr
Central age model (CAM) \cr
Finite mixture age model (FMM) \cr
Minimum age model with three parameters (MAM3) \cr
Minimum age model with four parameters (MAM4) \cr\cr
Function "SGCED" fits OSL dose-response curve and calculating equivalent dose using standardised growth curves (SGC) method, assiociated standard errors of equivlalent dose can also be estimated using either simple transformation or monte Carlo method.
}
\details{
\tabular{ll}{
Package: \tab RadialPlotter\cr
Type: \tab Package\cr
Version: \tab 1.7\cr
Date: \tab 2013-06-20\cr
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
\keyword{SGC method}
