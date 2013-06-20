\name{SGCED}
\alias{SGCED}
\alias{SGCED.default}
\title{
Fitting dose-response curve and calculating equivalent dose using standardised growth curves (SGC) method in optically stimulated luminescence (OSL) dating
}
\description{
Fitting OSL dose-response curve and calculating equivalent dose using standardised growth curves (SGC) method after Roberts HM and Duller GAT (2004), the error of equivalent dose is estimated after Duller GAT (2007). It is based on the fact that using information collected as part of standard SAR measurements, it is possible to construct reproducible patterns of growth from multiple-grain aliquots of a large number of samples (Roberts HM and Duller GAT, 2004).
}
\usage{
SGCED(CurveData, 
      Ltx, 
      model = "exp", 
      iniPars = NULL, 
      ErrorMethod = "mc", 
      MCiter = 10000, 
      plot = TRUE, 
      dashline = TRUE)
}
\arguments{
  \item{CurveData}{\link{vector}(\bold{required}): two columns, doses (dose1, dose2, ect) and standardized signals (Lx1/Tx1, Lx2/Tx2, ect) for determining the dose-response curve 
}
  \item{Ltx}{\link{vector}(\bold{required}): either one (standardized signal) or two columns (both standardized signal and its standard error), from which equivalent dose can be estimated by interpolating 
}
  \item{model}{\link{character}(with default): a model("linear", "exponential" or "linear+exponential") used for fitting the dose-response curve, see details
}
  \item{iniPars}{\link{numeric}(with default): initial parameters that will be used in dose-response curve fitting
}
  \item{ErrorMethod}{\link{character}(with default): method ("sp" or "mc") for estimating the standard error of equivalent dose if both standardized signal and its standard error are provided, see details
}
  \item{MCiter}{\link{numeric}(with default): maximum iterative number of the Monte Carlo method used for error assessing
}
  \item{plot}{\link{logical}(with default): whether drawing the dose-response curve or not
}
  \item{dashline}{\link{logical}(with default): whether add dash lines to the coordinate of estimated equivalent or not
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
  \item{errormessage}{error message generated during the calculation}
  \item{LMpars}{characteristic parameters and standard errors of the dose-response curve obtained by Levenberg-Marquardt method}
  \item{residual}{square of the total residual error}
  \item{fit.value}{observations .VS. fitted values}
  \item{Dose}{estimated equivalent dose (standard error)}
}
\references{
Roberts, H.M. and Duller, G.A.T., 2004. Standardised growth curves for optical dating of sediment using multiple-grain aliquots. Radiation Measurements 38(2), pp. 241-252.

Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from single aliquot regenerative dose measurements. Ancient TL 25(1), pp. 15-24.

Roberts, H.M., Durcan, J.A., Duller, G.A.T., 2009. Exploring procedures for the rapid assessment of optically stimulated luminescence range-finder ages. Radiation Measurements, 44(5-6), pp. 582-587.

Jorge More, Burton Garbow, Kenneth Hillstrom, User Guide for MINPACK-1, Technical Report ANL-80-74, Argonne National Laboratory, 1980.

\bold{Further reading} 

Long, H., Lai, Z.P., Fan, Q.S., Sun, Y.J., Liu, X.J., 2010. Applicability of a quartz OSL standardised growth curve for De determination up to
400 Gy for lacustrine sediments from the Qaidam Basin of the Qinghai-Tibetan Plateau. Quaternary Geochronology 5(2-3), pp. 212-217.

Yang, L.H., Lai, Z.P., Long, H., Zhang, J.R., 2011. Construction of a quartz OSL standardised growth curve (SGC) for aeolian samples from the Horqin dunefield in northeastern China. Geochronometria 38(4), pp. 391-396.
}
\author{
Peng Jun, email<pengjun10@mails.ucas.ac.cn>
}

\examples{
 #
 # Checking function SGC using data from Duller(2007)
 CurveData<-cbind(c(0,0.5,1,1.5,2.01),c(0.015,0.642,1.347,1.79,2.281),
                  c(0,10.2,20.5,30.7,41),c(0.011,2.395,4.078,5.295,6.209),
                  c(0,0.1,0.21,0.41,0.82),c(0.022,0.269,0.496,0.982,1.888),
                  c(0,0.31,0.62,1.24,2.47),c(0.09,0.521,0.953,1.655,3.364) )
  Ltx<-rbind(c(0.917,0.066),
             c(5.081,0.072),
             c(0.13,0.003),
             c(1.026,0.145) )
  Res<-matrix(nrow=4,ncol=3)
  for (i in 1:4)  {
    model<-ifelse(any(i==c(1,2)),"exp","line")    
    Simple<-SGCED(CurveData[,(2*i-1):(2*i),drop=FALSE],Ltx[i,,drop=FALSE],
                  model=model,ErrorMethod="sp",plot=FALSE)
    MonteCarlo<-SGCED(CurveData[,(2*i-1):(2*i),drop=FALSE],Ltx[i,,drop=FALSE],
                      model=model,ErrorMethod="mc",plot=FALSE)
    Res[i,1:2]<-Simple$Dose
    Res[i,3]<-MonteCarlo$Dose[1,2]
  } # end for
  colnames(Res)<-c("ED","Simp.Error","MC.Error")
  rownames(Res)<-c("Exponential","Exponential","linear","linear")
  print(Res)
 #
 #
 # that's my own data
 x<-cbind(c(0,18,36,54,72),c(0.03,1.49,2.51,3.32,4))
 Ltx<-cbind(c(3.117,3.321,2.925,2.998,2.913,2.825,3.007,2.972,3.082,3.083,
              2.999,3.737,2.74,3.13,3.77,3.15,2.12,2.35,3.01,2.93,2.77,3.15,2.33,1.98),
            c(0.1233,0.111,0.187,0.188,0.131,0.185,0.0972,0.193,0.124,0.112,
              0.196,0.2377,0.135,0.177,0.231,0.222,0.32,0.221,0.192,0.111,0.244,
              0.1811,0.197,0.207))
 SGCED(x,Ltx,model="line+exp")
}
\keyword{package RadialPlotter}
\keyword{OSL dating}
\keyword{SGC method}
\keyword{equivalent dose}
