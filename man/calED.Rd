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
calED(CurveData, 
      Ltx, 
      model = "line", 
      iniPars = NULL, 
      ErrorMethod = "mc", 
      MCiter = 1000, 
      plot = TRUE)
}
\arguments{
  \item{CurveData}{\link{vector}(\bold{required}): three columns, doses (dose1, dose2, ect), standardized signals (Lx1/Tx1, Lx2/Tx2, ect) 
and errors of standardized signals
}
  \item{Ltx}{\link{vector}(\bold{required}): two values (standardized signal and its standard error), from which equivalent dose can be estimated by interpolating
}
  \item{model}{\link{character}(with default): a model("linear", "exponential" or "linear+exponential") used for fitting the dose-response curve, see details
}
  \item{iniPars}{\link{numeric}(with default): initial parameters that will be used in dose-response curve fitting
}
  \item{ErrorMethod}{\link{character}(with default): method ("sp" or "mc") for estimating the standard error of equivalent dose, see details
}
  \item{MCiter}{\link{numeric}(with default): maximum iterative number if using Monte Carlo method for error assessing
}
  \item{plot}{\link{logical}(with default): whether drawing the dose-response curve or not
}
}
\details{
This function can used for both dose-response curve fitting and equivalent dose calculation, three fitting methods are available: \cr
1) "line": linear model, y=a*x+b; \cr
2) "exp": exponential model, y=a*(1-exp(-b*x))+c; \cr
3) "line+exp": exponential plus linear model, y=a*(1-exp(-b*x))+c*x+d \cr\cr
Standard error of equivalent dose can be accessed using two method: \cr
a) "sp": simple transformation of s(Lx/Tx) to s(ED); \cr
b) "mc": Monte Carlo method \cr\cr
Curve fitting is performed using the Levenberg-Marquardt algorithm (minpack: original FORTRAN 77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. FORTRAN 90 version by John Burkardt). Interpolation is performed using a combination of golden section search and successive parabolic interpolation (the same as function \link{optimize} in package "stats").
}
\value{
Returns a list that contains following elements:
  \item{mcED}{equivalent doses obtained by Monte Carlo simulation}
  \item{LMpars}{characteristic parameters and standard errors of the dose-response curve obtained by Levenberg-Marquardt method}
  \item{residual}{square of the total residual error}
  \item{fit.value}{observations .VS. fitted values}
  \item{ED}{estimated equivalent dose (standard error)}
}
\references{
Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from single aliquot regenerative dose measurements. 
Ancient TL 25(1), pp. 15-24.

Duller, G., 2007. Analyst. pp. 27-28.
}
\author{
Peng Jun, email<pengjun10@mails.ucas.ac.cn>
}
\note{
Interpolation is performed within the maximum re-dose, so the value of Ltx that used to calculate a equivalent dose must not larger than the maximum standardized signals that used to build a dose-response curve. 
}
\seealso{
  \link{sgcED}
}
\examples{
  #
  # load package 'RadialPlotter'
  library(RadialPlotter)
  #
  # test using data from Duller(2007)
  require(graphics)
  par(mfrow=c(2,2))
  CurveData<-cbind(c(0.0,0.5,1.0,1.5,2.01),c(0.015,0.642,1.347,1.79,2.281),
                   c(0.035,0.056,0.095,0.11,0.134),c(0.0,10.2,20.5,30.7,41.0),
                   c(0.011,2.395,4.078,5.295,6.209),c(0.00,0.034,0.059,0.076,0.089),
                   c(0.0,0.1,0.21,0.41,0.82),c(0.022,0.269,0.496,0.982,1.888),
                   c(0.002,0.005,0.008,0.016,0.029),c(0.0,0.31,0.62,1.24,2.47),
                   c(0.09,0.521,0.953,1.655,3.364),c(0.098,0.114,0.180,0.242,0.495))
  Ltx<-rbind(c(0.917,0.066),c(5.081,0.072),
             c(0.13,0.003),c(1.026,0.145) )
  Res<-matrix(nrow=4,ncol=3)
  for (i in 1:4)  {
    model<-ifelse(any(i==c(1,2)),"exp","line")    
    Simple<-calED(CurveData[,(3*i-2):(3*i),drop=FALSE],Ltx[i,],
                  model=model,ErrorMethod="sp",plot=FALSE)
    MonteCarlo<-calED(CurveData[,(3*i-2):(3*i),drop=FALSE],Ltx[i,],
                      model=model,ErrorMethod="mc")
    Res[i,1:2]<-Simple$ED
    Res[i,3]<-MonteCarlo$ED[2]
  } # end for
  colnames(Res)<-c("ED","Simp.Error","MC.Error")
  rownames(Res)<-c("Exponential","Exponential","linear","linear")
  print(Res)
  #
  # test using simulated data
  set.seed(123456789)
  ReDose<-rep(0:5,times=10)
  Ltx<-rep(c(0.015,0.642,1.347,1.79,2.281,2.53),10)+
       rnorm(60,mean=0,sd=c(0.005,0.01,0.05,0.1,0.1,0.1))
  sLtx<-rnorm(60,mean=0.2,sd=0.05)
  CurveData<-cbind(ReDose,Ltx,sLtx)
  calED(CurveData,c(2.136,0.195),model="exp",plot=FALSE)
}
\keyword{package RadialPlotter}
\keyword{OSL dating}
\keyword{equivalent dose}
\keyword{SAR method}
