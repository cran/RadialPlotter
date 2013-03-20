\name{RadialPlotter}
\alias{RadialPlotter}
\title{
 Estimating parameters of Galbraith's finite mixture age model and central age model and drawing radial plot  
}
\description{
Depending on the specified number of components, this function carrys out finite mixture age model analysis dynamically.
}
\usage{
RadialPlotter(EDdata, 
              ncomp = 0, 
              addsigma = 0, 
              maxiter = 500, 
              maxcomp = 9, 
              eps = .Machine$double.eps^0.5, 
              plot=TRUE,
              zscale = NULL, 
              samplename = NULL)
}
\arguments{
  \item{EDdata}{\link{data.frame}(\bold{required}): equivalent doses and associated errors (two columns), all the equivalent doses must be 
larger than 0
}
  \item{ncomp}{\link{numeric}(\bold{with default}): the number of components used to fit, 0 means analyzing automatical using finite mxiture age model, 1 for fitting the central age model, any other positive integer will fitting finite mixture age model with the specified number of components
}
  \item{addsigma}{\link{numeric}(\bold{with default}): the added spread to the relative error of Equivalent Dose (log-scale)
}
  \item{maxiter}{\link{numeric}(\bold{with default}): the maximum iterative number allowed for cnetral age model and finite mixture age model
}
  \item{maxcomp}{\link{numeric}(\bold{with default}): the maximum number of components allowed for fitting finite mixture age model, which has
no up limitation, but it will be time consuming and make no sense if it is too large
}
  \item{eps}{\link{numeric}(\bold{with default}): the maximum tolerance for stopping the iterative process, iterative process will stop if updated parameters satisfying ||pars-newpars||<=eps
}
  \item{zscale}{\link{numeric}(\bold{optional}): option for modifying the z-scale of radial plot manually. Example: \code{zscale = seq(min(EDdata),max(EDdata), by = 3L)}
}
  \item{plot}{\link{logical}(\bold{with default}): if drawing the radial plot using estimated parameters
}
  \item{samplename}{\link{character}(\bold{optional}): the name of the sample used for fitting
}
}
\details{
This function fits finite mixture age model with a given component number (by specifying a number of components that not smaller than 2) or do the fitting automatically (if ncomp=0), or analyzes the central age model if the number of components is specified to be 1, also Notes that all these calculations are performed using log-scale, which means that any minus equivalent dose is not allowed. 
}
\value{
Returns a radial plot, a terminal output containing estimated results, accompanying with a list that contains following elements:

  \item{errorflag}{0 indicates a successful work, 1 implicates that the covariance matrix is singular or the diagonal position of the inversed covariance matrix have at least one minus value, which typically indicates an incorrect addsigma value or an inappropriate number of components
}
  \item{pars}{the estimated parameters and standard errors}
  \item{BIC}{the BIC value}
  \item{maxlik}{the maximum logged likelihood value}
}
\references{
Galbraith, R.F., 1988. Graphical Display of Estimates Having Differing Standard Errors. Technometrics, 30 (3), pp. 271-281.

Galbraith, R.F., 1990. The radial plot: Graphical assessment of spread in ages. International Journal of Radiation Applications and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 17 (3), pp. 207-214.

Galbraith, R.F., Green, P., 1990. Estimating the component ages in a finite mixture. International Journal of Radiation Applications and Instrumentation. Part D. Nuclear Tracks and Radiation Measurements, 17 (3), pp. 197-206.

Galbraith, R.F., Laslett, G.M., 1993. Statistical models for mixed fission track ages. Nuclear Tracks And Radiation Measurements, 21 (4), pp. 459-470.

Galbraith, R.F., 1994. Some Applications of Radial Plots. Journal of the American Statistical Association, 89 (428), pp. 1232-1242.

Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley, J.M., 1999. Optical dating of single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design and statistical models. Archaeometry, 41, pp. 339-364.

Galbraith, R.F., 2005. Statistics for Fission Track Analysis, Chapman & Hall/CRC, Boca Raton.

Galbraith, R.F., 2010. On plotting OSL equivalent doses. Ancient TL, 28 (1), pp. 1-10.

Galbraith, R.F., Roberts, R.G., 2012. Statistical aspects of equivalent dose and error calculation and display in OSL dating: An overview and some recommendations. Quaternary Geochronology,11, pp.1-27.

Sebastian, K., Christoph, S., Margret, C., F., Michael, D., Manfred, F., Markus, F., 2012. Introducing an R package for luminescence dating analysis. Ancient TL, 30, pp. 1-8.

\bold{Further reading}

Vermeesch, P., 2009. RadialPlotter: a Java application for fission track, luminescence and other radial plots, Radiation Measurements, 
44 (4), pp. 409-410.

}

\author{Peng Jun}

\note{The original S code for drawing radial plot was written by Rex Gablraith and was transformed to R code by Sebastian Kreutzer. The code for drawing radial plot in this package is modified from package 'Luminescence' by Sebastian Kreutzer et al (2012).
}

\examples{
#
# loading the equivalent dose data
  data(EDdata)
#
# specifying the number of component to be 1 to do central age model analysis
  RadialPlotter(EDdata,ncomp=1,plot=FALSE)
#
# specifying the number of component to be 2 to fitting the finite mxiture age
# model with  a number of components of 2 (k=2)
  RadialPlotter(EDdata,ncomp=2)
#
# analyzing equivalent dose data by using the default option to pick
# out the appropriate number of components then fitting automatically
  RadialPlotter(EDdata,maxcomp=13,samplename='Tengger')
#
# setting zscale manually
  zscale = seq(16,max(EDdata), by = 3L)
  RadialPlotter(EDdata,zscale=zscale)
}

\keyword{RadialPlotter}

