\name{decomp}
\alias{decomp}
\alias{decomp.default}
\title{
Decomposition decay curve in optically stimulated luminescence (OSL) dating
}
\description{
Function "decomp" decomposes the OSL decay curve to a specified number of components using a combination of differential evolution and
Levenberg-Marquardt algorithm.
}
\usage{
decomp(SigData, 
       ncomp = 3, 
       LEDpower = 60, 
       LEDwavelength = 470, 
       plot = TRUE, 
       xyLim = c(5, 500), 
       typ = "CW", 
       control.args = list())
}

\arguments{
  \item{SigData}{\link{data.frame}(\bold{required}): two columns, time (t1,t2,ect) and correspond signal (s1,s2,ect) to be decomposed
}
  \item{ncomp}{\link{numeric}(with default): the number of components to be decomposed, which must be in a range of 1 to 7
}
  \item{LEDpower}{\link{numeric}(with default): LED power (mW/cm2) used for intensity ramping during the measurement

}
  \item{LEDwavelength}{\link{numeric}(with default): LED wavelength (nm) used for stimulation during the measurement
}
  \item{plot}{\link{logical}(with default): whether plot the results or not
}
  \item{xyLim}{\link{vector}(with default): up limitation of x and y axis used for plotting
}
  \item{typ}{\link{character}(with default): type of the OSL decay curve, only type "CW" is avialable currently
}
  \item{control.args}{\link{list}(with default): further arguments used in differential evolution algorithm, see details
}
}
\details{
This function decomposes OSL decay curve to a specified number of components using a combination of differential evolution and
Levenberg-Marquardt algorithm. Fitting model: \cr
I(t)=a1*exp(-b1*t)+a2*exp(-b2*t)+...+ak*exp(-bk*t) \cr\cr
Initial parameters are estimated using the differential evolution method suggested by Bluszcz A and Adamiec G (2006), then Levenberg-Marquardt algorithm (minpack: original FORTRAN 77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. FORTRAN 90 version by John Burkardt) will be performed to obtain a final optimization. In some radical cases, differential evolution will fail in parameter initializing, then a serise of simple trails will be carryed out using Levenberg-Marquardt algorithm. The arguments (control.args) used for differential evolution are:\cr\cr
factor: used for control the number of parents (NP), np=factor*ncomp, 5 to 10 times the number of components is preferable, default factor=10; \cr\cr
f: weighting factor, ranging from 0 to 1.2, default f=0.5; \cr\cr
cr: crossover constant, ranging from 0 to 1, default cr=0.99; \cr\cr
maxiter: the maximum iterative number allowed, default maxiter=1e3; \cr\cr
tol: a tolerance for stopping the iteration, if magnitude of all standard deviations of parameters within the whole population are smaller than
tol, the iteration will be terminated. default tol=0.1
}
\value{
Returns a list containing following elements:\cr
  \item{Comp.Signal}{fitted values for observations as well as for each component}
  \item{pars}{estimated parameters and associated standard errors}
  \item{value}{square of the total residual error}
  \item{errorflag}{error message generated during the calculation, if parameters' standard errors can be estimated, errorflag=0, else 1}
}
\references{
Bluszcz, A., 1996. Exponential function fitting to TL growth data and similar applications. Geochronometria, 13, pp.135-141.

Bluszcz, A., Adamiec, G., 2006. Application of differential evolution to fitting OSL decay curves. Radiation Measurements, 41 (7-8), pp.886-891.
 
Differential evolution, http://en.wikipedia.org/wiki/Differential_evolution
  
Jain, M., Murray, A.S., Boetter-Jensen, L., 2003. Characterisation of blue-light stimulated luminescence components in different quartz samples: implications for dose measurement. Radiation Measurements, 37 (4-5), pp.441-449.

Jorge, M., Burton, G., Kenneth, H., 1980. User Guide for MINPACK-1, Technical Report ANL-80-74, Argonne National Laboratory. 
}
\author{
Peng Jun, email<pengjun10@mails.ucas.ac.cn>
}
\note{Both LEDpower and LEDwavelength are used for determining the absolute photoionisation cross section, hence to achieve desirable results
for a sample of interest, making sure correct LEDpower and LEDwavelength have been employed.\cr\cr
In most cases, parameters can be initialized by differential evolution method, but sometimes Levenberg-Marquardt algorithm fails to obtain
a final optimization (that's exactly the truth as the increasing of the number of components, which might indicate overcharging of the model), if this happened, parameters initialized with differential evolution method will be regarded as the final estimation.\cr\cr
I thank Andrzej Bluszcz, who helps me a lot during the writting of this function.
}

\examples{
 #
 # loading signal data
 data(CW.Signal.data)
 #
 # decompose OSL signal to 3 components
 decomp(CW.Signal.data,ncomp=3)
 #
 # reset factor, maxiter and decompose signal to 5 components
 decomp(CW.Signal.data,ncomp=5,control.args=list(factor=9,maxiter=500))
}

\keyword{package RadialPlotter}
\keyword{OSL dating}
\keyword{decay curve decomposition}


