\name{calED}
\alias{calED}
\alias{calED.default}
\title{
Fitting dose-response curve and calculating equivalent dose in optically stimulated luminescence
}
\description{
Fitting OSL dose-response curve and calculating equivalent dose, the error of equivalent dose is estimated after Duller GAT (2007).
}
\usage{
calED(Curvedata, Ltx, 
      model = c("line","exp","line+exp"), origin = FALSE,
      nstart = 100, upb = 1, ErrorMethod = c("mc", "sp"),
      nsim = 1000, plot = TRUE, samplename = NULL)
}
\arguments{
  \item{Curvedata}{\link{data.frame}(\bold{required}): three columns, regenerative doses (dose1, dose2, ect), standardized signals (Lx1/Tx1, Lx2/Tx2, ect) 
and errors of standardized signals
}
  \item{Ltx}{\link{vector}(\bold{required}): two values (standardized natural signal and its standard error), from which equivalent dose and its standard error can be estimated with interpolation
}
  \item{model}{\link{character}(with default): a model ("line", "exp" or "line+exp") used for fitting the dose-response curve, default model="line". See details
}
  \item{origin}{\link{logical}(optional): whether force the fitting to pass the origin (x=0,y=0) or not
}
  \item{nstart}{\link{numeric}(with default): maximum number of attempts that used to initialize parameters in curve fitting
}
  \item{upb}{\link{numeric}(with default): upper boundary of b value, initial b values will be generated uniformly from (0, upb), see details
}
  \item{ErrorMethod}{\link{character}(with default): method ("sp" or "mc") for estimating the standard error of equivalent dose. Default ErrorMethod="mc". See details for more information
}
  \item{nsim}{\link{numeric}(with default): maximum simulative number allowed if using Monte Carlo method for error assessing
}
  \item{plot}{\link{logical}(with default): whether drawing the fitted dose-response curve and calculated equivalent dose or not
}
  \item{samplename}{\link{character}(with default): name of the sample
}
}
\details{
This function can used for both dose-response curve fitting and equivalent dose calculating, three fitting methods are available: \cr
1) "line": linear model, y=a*x+b \bold{or} y=a*x if origin=TRUE;\cr
2) "exp": exponential model, y=a*(1-exp(-b*x))+c \bold{or} y=a*(1-exp(-b*x)) if origin=TRUE;\cr
3) "line+exp": exponential plus linear model, y=a*(1-exp(-b*x))+c*x+d \bold{or} y=a*(1-exp(-b*x))+c*x if origin=TRUE.\cr\cr
Standard error of equivalent dose can be accessed using two method (Duller GAT, 2007): \cr
a) "sp": simple transformation of s(Lx/Tx) to s(ED); \cr
b) "mc": Monte Carlo method \cr\cr
Curve fitting is performed using the Levenberg-Marquardt algorithm (minpack: original FORTRAN 77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. FORTRAN 90 version by John Burkardt). Interpolation is performed using a combination of golden section search and successive parabolic interpolation (the same as function \link{optimize} in package "stats").\cr\cr
As for fitting a linear model, initial parameters will not be a problem, but if the model to be fitted is either exponential or exponential plus linear, then a series of b values will be initialized randomly from a uniform distribution in the range of 0 to upb, hence other parameters (a or a, c or a, c, d) can be estimated with Linear Algebra method. User is advised to keep the default setting "plot" to be TRUE if possible in order to see whether the fitted dose-response curve is of the correct form.
}
\value{
Return an invisible list that contains following elements:
  \item{mcED}{equivalent doses obtained by Monte Carlo simulation, available only if ErrorMethod="mc", else it will be NULL}
  \item{LMpars}{characteristic parameters and standard errors of the dose-response curve estimated with Levenberg-Marquardt method}
  \item{residual}{sum of square of residual errors}
  \item{fit.value}{observations .VS. fitted values}
  \item{ED}{estimated equivalent dose (and corresponded standard error)}
}
\references{
Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from single aliquot regenerative dose measurements. 
Ancient TL 25(1), pp. 15-24.

Duller, G., 2007. Analyst. pp. 27-28.
}
\note{
Interpolation is performed within the maximum regenerative dose, so the value of Lt/Tx that used to calculate an equivalent dose must not larger than the maximum standardized signals that used to build the dose-response curve. 
}
\seealso{
  \link{sgcED}; \link{optimize}
}
\examples{
  # load package "RadialPlotter"
  library(RadialPlotter)
  # example 1: test using data from Duller (2007)
  x<-data.frame(c(0.0,0.5,1.0,1.5,2.01),c(0.015,0.642,1.347,1.79,2.281),
                c(0.035,0.056,0.095,0.11,0.134),c(0.0,10.2,20.5,30.7,41.0),
                c(0.011,2.395,4.078,5.295,6.209),c(0.00,0.034,0.059,0.076,0.089),
                c(0.0,0.1,0.21,0.41,0.82),c(0.022,0.269,0.496,0.982,1.888),
                c(0.002,0.005,0.008,0.016,0.029),c(0.0,0.31,0.62,1.24,2.47),
                c(0.09,0.521,0.953,1.655,3.364),c(0.098,0.114,0.180,0.242,0.495))
  Ltx<-rbind(c(0.917,0.066),c(5.081,0.072),
             c(0.13,0.003),c(1.026,0.145))
  Res<-matrix(nrow=4,ncol=3)
  for (i in 1:4)  {
    model<-ifelse(any(i==c(1,2)),"exp","line")    
    Simple<-calED(x[,(3*i-2):(3*i),drop=FALSE],Ltx[i,],
                  model=model,ErrorMethod="sp",plot=FALSE)
    MonteCarlo<-calED(x[,(3*i-2):(3*i),drop=FALSE],Ltx[i,],
                      model=model,ErrorMethod="mc",plot=FALSE)
    Res[i,1:2]<-Simple$ED
    Res[i,3]<-MonteCarlo$ED[2]
  } # end for
  colnames(Res)<-c("ED","Simp.Error","MC.Error")
  rownames(Res)<-c("Exponential","Exponential","linear","linear")
  print(Res)
  # example 2
  Curvedata<-data.frame(c(0,18,36,54,72,0,18),                        # ReDose
                        c(0.0272,1.494,2.51,3.325,4.0,0.0272,1.453),  # Lx/Tx
                        c(0.005,0.067,0.097,0.109,0.201,0.005,0.053)) # s(Lx/Tx)
  Res<-calED(Curvedata,c(3.117,0.1),model="line+exp",plot=TRUE,origin=TRUE)
  Res
}
\keyword{package RadialPlotter}
\keyword{OSL dating}
\keyword{equivalent dose}
\keyword{SAR method}
