\name{sgcED}
\alias{sgcED}
\alias{sgcED.default}
\title{
Fitting dose-response curve and calculating equivalent dose using standardised growth curves (SGC) method
}
\description{
Its a wrapped version of function "calED', fitting OSL dose-response curve and calculating equivalent dose values using standardised growth curves (SGC) method after Roberts HM and Duller GAT (2004), the error of equivalent doses are estimated after Duller GAT (2007). It is based on the fact that using information collected as part of standard SAR measurements, it is possible to construct reproducible patterns of growth from multiple-grain aliquots of a large number of samples (Roberts HM and Duller GAT, 2004).
}
\usage{
sgcED(Curvedata, 
      Ltx, 
      model = c("line","exp","line+exp"),
      nstart = 100,
      upb = 5,
      ErrorMethod = c("mc","sp"),
      MCiter = 1000, 
      plot = TRUE,
      samplename = NULL,
      outfile = NULL)
}
\arguments{
  \item{Curvedata}{\link{data.frame}(\bold{required}): three columns, the same as that in function \link{calED}
}
  \item{Ltx}{\link{data.frame}(\bold{required}): two columns (standardized signals and its standard errors), from which equivalent doses can be estimated by interpolating 
}
  \item{model}{\link{character}(with default): a model("linear", "exponential" or "linear+exponential") used for fitting the dose-response curve, the same as that in function \link{calED}
}
  \item{nstart}{\link{numeric}(with default): maximum number of attempts that used to initialize parameters in curve fitting, the same as that in function \link{calED}
}
  \item{upb}{\link{numeric}(with default): upper boundary of b value, initial b values will be generated uniformly from (0, upb), the same as that in function \link{calED}
}
  \item{ErrorMethod}{\link{character}(with default): method ("sp" or "mc") for estimating the standard errors of equivalent doses, 
see \link{calED} for more information
}
  \item{MCiter}{\link{numeric}(with default): maximum iterative number if using Monte Carlo method for error assessing
}
  \item{plot}{\link{logical}(with default): whether drawing a plot or not
}
  \item{samplename}{\link{character}(optional): name of the sample
}
  \item{outfile}{\link{character}(optional): if specified, calculated equivalent dose values will be written to a file of name
outfile in .csv format}
}
\details{
see function \link{calED} for details.
}
\value{
  Return a visible list that contains following elements:
  \item{LMpars}{characteristic parameters and standard errors of the dose-response curve obtained by Levenberg-Marquardt method}
  \item{residual}{square of the total residual error}
  \item{fit.value}{observations .VS. fitted values}
  \item{ED}{estimated equivalent doses (standard errors)}
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
\seealso{
  \link{calED}
}
\examples{
 # example 1
 x<-data.frame(c(0,18,36,54,72),                  # Re-Dose of dose-response curve
               c(0.03,1.49,2.51,3.32,4),          # Lx/Tx of dose-response curve
               c(0.001,0.011,0.016,0.02,0.026))   # sLtx
 #
 Ltx<-data.frame(c(3.117,3.321,2.925,2.998,2.913,2.825,3.007,         # Lt/Tx 
                   2.972,3.082,3.083,2.999,3.737,2.74,3.13,3.77,
                   3.15,2.12,2.35,3.01,2.93,2.77,3.15,2.33,1.98),    
                 c(0.1233,0.111,0.187,0.188,0.131,0.185,0.0972,       # sLtx 
                   0.193,0.124,0.112, 0.196,0.2377,0.135,0.177,0.231,
                   0.222,0.32,0.221,0.192,0.111,0.244,0.1811,0.197,0.207)) 
 sgcED(x,Ltx,model="line+exp")
 # example 2
 ltx1<-c(0.032,1.61,2.55,3.21,3.87,0.031,1.55)
 ltx2<-c(0.025,1.44,2.47,3.35,4.17,0.033,1.47)
 ltx3<-c(0.027,1.51,2.68,3.52,4.41,0.021,1.39)
 ltx<-cbind(ltx1,ltx2,ltx3)
 ltx<-cbind(apply(ltx,MARGIN=1,mean),apply(ltx,MARGIN=1,sd))
 redose<-c(0,18,36,54,72,0,18)
 Curvedata<-data.frame(redose,ltx)
 Ltx<-data.frame(seq(0.5,4.0,by=0.5),rep(0.1,8))
 sgcED(Curvedata,Ltx,model="line+exp")
}
\keyword{package RadialPlotter}
\keyword{OSL dating}
\keyword{equivalent dose}
\keyword{SGC method}
