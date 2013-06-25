\name{sgcED}
\alias{sgcED}
\alias{sgcED.default}
\title{
Fitting dose-response curve and calculating equivalent dose using standardised growth curves (SGC) method
}
\description{
Its a wrapped version of function "calED', fitting OSL dose-response curve and calculating equivalent dose using standardised growth curves (SGC) method after Roberts HM and Duller GAT (2004), the error of equivalent doses are estimated after Duller GAT (2007). It is based on the fact that using information collected as part of standard SAR measurements, it is possible to construct reproducible patterns of growth from multiple-grain aliquots of a large number of samples (Roberts HM and Duller GAT, 2004).
}
\usage{
sgcED(CurveData, 
      Ltx, 
      model = "line", 
      iniPars = NULL, 
      ErrorMethod = "mc", 
      MCiter = 1000, 
      plot = TRUE)
}
\arguments{
  \item{CurveData}{\link{vector}(\bold{required}): three columns, the same as in function \link{calED}
}
  \item{Ltx}{\link{vector}(\bold{required}): two columns (standardized signals and its standard errors), from which equivalent doses can be estimated by interpolating 
}
  \item{model}{\link{character}(with default): a model("linear", "exponential" or "linear+exponential") used for fitting the dose-response curve, the same as in function \link{calED}
}
  \item{iniPars}{\link{numeric}(with default): initial parameters that will be used in dose-response curve fitting
}
  \item{ErrorMethod}{\link{character}(with default): method ("sp" or "mc") for estimating the standard errors of equivalent doses, 
see \link{calED}
}
  \item{MCiter}{\link{numeric}(with default): maximum iterative number if using Monte Carlo method for error assessing
}
  \item{plot}{\link{logical}(with default): whether drawing a plot or not
}
}
\details{
see function \link{calED} for details.
}
\value{
  Returns a list that contains following elements:
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
 # CurveData
 x<-cbind(c(0,18,36,54,72),                # Re-Dose of dose-response curve
          c(0.03,1.49,2.51,3.32,4),        # Lx/Tx of dose-response curve
          c(0.001,0.011,0.016,0.02,0.026)) # sLtx
 #
 # Ltx used for calculate equivalent doses and corresponded standard errors
 Ltx<-cbind(c(3.117,3.321,2.925,2.998,2.913,2.825,3.007,2.972,3.082,3.083,2.999,
              3.737,2.74,3.13,3.77,3.15,2.12,2.35,3.01,2.93,2.77,3.15,2.33,1.98),    # Lt/Tx 
            c(0.1233,0.111,0.187,0.188,0.131,0.185,0.0972,0.193,0.124,0.112, 0.196,0.2377,
              0.135,0.177,0.231,0.222,0.32,0.221,0.192,0.111,0.244,0.1811,0.197,0.207)) # sLtx 
 sgcED(x,Ltx,model="line+exp")
}
\keyword{package RadialPlotter}
\keyword{OSL dating}
\keyword{equivalent dose}
\keyword{SGC method}
