###     The R package RadialPlotter
###     R Package for OSL equivalent dose analysis and radial plot drawing.
###     Copyright (C) 2013 Peng Jun 
###
###     Written by Peng Jun.
###
###     This file the R scripts (all R functions) for package RadialPlotter.
###     URL: http://CRAN.R-project.org/package=RadialPlotter
### 
###     The RadialPlotter package is free software: you can redistribute it  
###     or modify it under the terms of the GNU General Public License as 
###     published by the Free Software Foundation, either version 3 of the 
###     License, or any later version.
###
###     Package RadialPlotter is distributed in the hope that it will be 
###     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
###     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###     GNU General Public License for more details.
###
###     You should have received a copy of the GNU General Public License
###     along with RadialPlotter.  If not, see <http://www.gnu.org/licenses/>.
###
###     If you wish to report bugs please contact the author:
###
###     Peng Jun
###     pengjun10@mails.ucas.ac.cn
###     University of Chinese Academy of Sciences
############################################### FUNCTION RadialPlotter ###############################################
### ******************************************************************************************************************
#### Function RadialPlotter is used to perform Galbraith's
#### statistical age models analysis. Models that can be 
###  analyzed include: CAM, FMM and MAM.
###
###     Author: Peng Jun, 2013.04.20, revised in 2013.07.26, revised in 2013.09.19
###
### References: Galbraith, R.F., 1988. Graphical Display of Estimates Having Differing Standard 
###             Errors. Technometrics, 30 (3), pp. 271-281.
###
###             Galbraith, R.F., Roberts, R.G., 2012. Statistical aspects of equivalent dose and error calculation and  
###             display in OSL dating: An overview and some recommendations. Quaternary Geochronology, 11, pp. 1-27.
### ******************************************************************************************************************
RadialPlotter<-
function(EDdata,ncomp=0,addsigma=0,maxiter=500,
         maxcomp=9,algorithm=c("lbfgsb","port"),
         eps=.Machine$double.eps^0.5,plot=TRUE,
         pcolor="blue",psize=1.5,kratio=0.3,
         zscale=NULL,samplename=NULL)  {
  UseMethod("RadialPlotter")
}
### set default method for function RadialPlotter
RadialPlotter.default<-
function(EDdata,ncomp=0,addsigma=0,maxiter=500,
         maxcomp=9,algorithm=c("lbfgsb","port"),
         eps=.Machine$double.eps^0.5,plot=TRUE,
         pcolor="blue",psize=1.5,kratio=0.3,
         zscale=NULL,samplename=NULL)  {
  ### stop if not 
  if(!is.data.frame(EDdata))    stop("Error: EDdata must be of type data.frame!")
  if(ncol(EDdata)!=2L)          stop("Error: EDdata must contain two columns!")
  if(!is.numeric(
     unlist(unclass(EDdata))))  stop("Error: elements in EDdata must be all of type numeric!")
  if(any(EDdata[,1L]<=0))       stop("Error: equivalent dose must be larger than zero!")
  if(any(EDdata[,2L]<=0))       stop("Error: std.error of equivalent dose must be larger than zero!")
  if(!is.numeric(maxcomp))      stop("Error: maxcomp must be of type numeric!")
  if(length(maxcomp)!=1L)       stop("Error: maxcomp must be an one-element vector!")
  if(abs(maxcomp-round(maxcomp))
     >=.Machine$double.eps^0.5) stop("Error: maxcomp must be a integer!")                                
  if(maxcomp<0)                 stop("Error: maxcomp must not be smaller than 0!")
  if(maxcomp>nrow(EDdata))      stop("Error: maxcomp must not exceed the length of EDdata!")
  if(!is.numeric(ncomp))        stop("Error: ncomp must be of type numeric!")
  if(length(ncomp)!=1L)         stop("Error: ncomp must be an one-element vector!")
  if(!ncomp %in% (-2:maxcomp))  stop("Error: ncomp must be an integer ranging from -2 to maxcomp!")
  if(!is.numeric(addsigma))     stop("Error: addsigma must be of type numeric!")
  if(length(addsigma)!=1L)      stop("Error: addsigma must be an one-element vector!")
  if(addsigma<0)                stop("Error: addsigma must not be smaller than 0!")
  if(!is.numeric(maxiter))      stop("Error: maxiter must be of type numeric!")
  if(length(maxiter)!=1L)       stop("Error: maxiter must be an one-element vector!")
  if(abs(maxiter-round(maxiter))
     >=.Machine$double.eps^0.5) stop("Error: maxiter must be an integer!")
  if(maxiter<10L)               stop("Error: maxiter is too small!")
  if(maxiter>10000L)            stop("Error: maxiter is too large!")
  if(!is.character(algorithm))  stop("Error: algorithm must be of type character!")
  if(length(algorithm)==1L) {
    if(!algorithm %in% c("lbfgsb","port"))      stop("Error: algorithm must be either 'lbfgsb' or 'port'!")
  } else if(length(algorithm)==2L) {
    if(!all(algorithm %in% c("lbfgsb","port"))) stop("Error: only algorithm 'lbfgsb' or 'port' is available!")
  } else {                                      stop("Error: algorithm must contain no more than two elements!")
  } # end if
  if(!is.numeric(eps))          stop("Error: eps must be of type numeric!")
  if(length(eps)!=1L)           stop("Error: eps must be an one-element vector!")
  if(eps<0)                     stop("Error: eps must be not smaller than 0!")
  if(!is.logical(plot))         stop("Error: plot must be either TRUE or FALSE!")
  if(length(plot)!=1L)          stop("Error: plot must be an one-element vector!")
  if(!is.character(pcolor))     stop("Error: pcolor must be of type character!")
  if(length(pcolor)!=1L)        stop("Error: pcolor must be an one-element vector!")
  if(!is.numeric(psize))        stop("Error: psize must be of type numeric!")
  if(length(psize)!=1L)         stop("Error: psize must be an one-element vector!")
  if(psize<0.1)                 stop("Error: psize is too small!")
  if(psize>5.0)                 stop("Error: psize is too large!")
  if(!is.numeric(kratio))       stop("Error: kratio must be of type numeric!")
  if(length(kratio)!=1L)        stop("Error: kratio must be an one-element vector!")
  if(abs(kratio)>3.0)           stop("Error: magnitude of kratio is too large!")
  if(!is.null(zscale) &&
     !is.numeric(zscale))       stop("Error: zscale must either be NULL or a numeric vector!")
  if(!is.null(samplename) &&
     !is.character(samplename)) stop("Error: samplename must either be NULL or of type character!")
  if(ncomp==0L && maxcomp<=1L)  stop("Error: maxcomp must non't be smaller than 2 if ncomp is 0!")  
  ### set n, ED, sED
  n<-nrow(EDdata) 
  ed<-drop(EDdata[,1L])
  error<-drop(EDdata[,2L])
  ### function for radial plot drawing,
  ### the code is revised from Rex Galbraith 
  ### *************************************************
  RadialPlot<-function(Data,Pars,zscale,samplename,kratio,pcolor,psize)  {
    z.i<-log(Data[,1L])
    se.i<-Data[,2L]/Data[,1L]
    ### if Pars=NULL, lines will not be plotted out
    if(!is.null(Pars))  {
      Pars<-log(Pars)
    } # end if
    plot.area_ratio<-kratio
    z.0<-sum(z.i/se.i^2)/sum(1/se.i^2)
    x.i<-1/se.i
    y.i<-(z.i-z.0)/se.i
    h<-plot.area_ratio*(max(x.i)-min(x.i))/(max(y.i)-min(y.i))
    r.0<-max(sqrt(x.i^2+y.i^2*h^2))
    if(is.null(zscale))  {
      mkest<-round(exp(pretty(z.i)), 0L)
      if(sd(z.i)/mean(z.i)<0.05)  {
        mkest<-round(exp(pretty(c(min(z.i)*0.98, max(z.i)*1.02))), 0L)
      } # end if
    } else  {
      mkest<-zscale
    } # end if
    mkr<-range(mkest)
    circ.x<-(0:200)/200
    circ.x<-mkr[1L]+circ.x*(mkr[2L]-mkr[1L])
    circ.z<-log(circ.x)
    zmkest<-log(mkest) 
    circ.x<-r.0/sqrt(1+h^2*(circ.z-z.0)^2)
    circ.y<-(circ.z-z.0)*circ.x
    xmk1<-r.0/sqrt(1+h^2*(zmkest-z.0)^2)
    xmk2<-xmk1*1.01
    xmk3<-1.01*1.04*r.0/sqrt(1+h^2*(zmkest-z.0)^2)
    ymk1<-(zmkest-z.0)*xmk1
    ymk2<-(zmkest-z.0)*xmk2
    ymk3<-(zmkest-z.0)*xmk3 
    may<-max(abs(y.i))
    yaxis.max<-if(may>12)   {
      may  } else if(may<3) {
      8    } else  {
      12   } # end if
    par(mfrow=c(1,1), oma=c(2,0,0,0), mar=c(5,4,4,4.5), xpd=TRUE, las=1, cex=1)
    plot(NA, NA, xlim=c(0,max(x.i)), ylim=c(-yaxis.max,yaxis.max), xlab="", ylab="",
       xaxt="n", yaxt="n", main=samplename, bty="n", typ="n", xaxs="i", yaxs="i")      
    lenpar<-length(Pars) 
    maxx<-which.max(circ.x)
    ### if Pars!=NULL, starting drawing lines out
    if(!is.null(Pars))  {
      for(i in 1:lenpar)  {
        if(Pars[i]-z.0>=0)  {
          Ci<-approx(circ.y[maxx:201], circ.x[maxx:201], ties="ordered",
          xout=(Pars[i]-z.0)*max(x.i), rule=2, method="constant", f=0.5)$y
        }  else  {
          Ci<-approx(circ.y[1:(maxx-1)], circ.x[1:(maxx-1)], ties="ordered",
          xout=(Pars[i]-z.0)*max(x.i), rule=2, method="constant",f=0.5)$y
        } # end if
       lines(c(0,Ci), c(0,(Pars[i]-z.0)*Ci), lty=1, col="black", lwd=1.5)
      } # end for
    } # end if
    lines(circ.x, circ.y, lwd=2)
    segments(xmk1, ymk1, xmk2, ymk2)
    text(xmk3, ymk3, round(mkest,2L), cex=1)
    text(max(x.i)*1.15,0,"De (Gy)", cex=1, srt=90)
    par(mfrow=c(1,1), oma=c(2,0,0,0.5), mar=c(5,4,4,4.5), xpd=TRUE, las=1, cex=1, new=TRUE)
    plot(NA, NA, xlim=c(0,max(x.i)), ylim=c(-yaxis.max, yaxis.max), yaxt="n", xaxt="n",
         xaxs="i", yaxs="i", ylab="Standardised Estimate", xlab="", type="n", bty="n")
    points(x.i, y.i, pch=21, col="black", cex=psize, bg=pcolor)
    axis(side=1, line=4, cex.axis=1, lwd=2)
    mtext(side=1, line=6, "Precision", cex=1)
    reticks.labels<-round(1/axTicks(side=1)*100,1L)
    reticks.values<-axTicks(side=1)
    axis(side=1, at=reticks.values[-1L], lwd=2, labels=reticks.labels[-1L], line=4, cex.axis=1, tck=0.02, padj=-4)
    mtext(side=1, line=1.5,"Relative Error (%)", cex=1)
    axis(side=2,at=c(-2,-1,0,1,2), lwd=2, labels=c("-2","","0","","2"), cex.axis=1)
    on.exit(par(oma=c(0,0,0,0),
                mar=c(5,4,4,2)+0.1))
  } ### end function RadialPlot
  ### ******************************************************
  tol<-.Machine$double.eps^0.3     
  ### function Rmam for port routines to do
  ### optimization for MAM3 or MAM4, added in 2013.07.25
  ### ******************************************************
  Rmam<-function(EDdata,ncomp,addsigma,maxiter) {
    ### -1 for 3-parameter, -2 for 4-parameter
    stopifnot(ncomp%in%c(-1L,-2L))
    ### add spread dispersion to relatvie 
    ### standard error of Equivalent Dose
    x<-sqrt((EDdata[,2L]/EDdata[,1L])^2+addsigma^2)
    ### change Equivalent Dose to log-scale
    y<-log(drop(EDdata[,1L]))
    ### function minfunc to be minimized
    ### ******************************
    minfunc<-function(p)  {
    ### here ncomp is a gloable variable
    if(ncomp==-1L) {
      ### for MAM3
      u0<-(p[2L]/(p[3L])^2+y/x^2)/(1/(p[3L])^2+1/x^2)
      sigma0<-1/(sqrt(1/(p[3L])^2+1/x^2)) 
      prop1<-p[1L]/sqrt(2*pi)/x*exp(-(y-p[2L])^2/2/x^2)
      prop2<-(1-p[1L])/sqrt(2*pi*((p[3L])^2+x^2))* 
             (1-pnorm((p[2L]-u0)/sigma0))/(1-pnorm((p[2L]-p[2L])/p[3L]))*
              exp(-(y-p[2L])^2/2/((p[3L])^2+x^2))
      return(-sum(log(prop1+prop2)))
    } else if(ncomp==-2L) {
      ### for MAM4
      u0<-(p[3L]/(p[4L])^2+y/x^2)/(1/(p[4L])^2+1/x^2)
      sigma0<-1/(sqrt(1/(p[4L])^2+1/x^2)) 
      prop1<-p[1L]/sqrt(2*pi)/x*exp(-(y-p[2L])^2/2/x^2)
      prop2<-(1-p[1L])/sqrt(2*pi*((p[4L])^2+x^2))* 
             (1-pnorm((p[2L]-u0)/sigma0))/(1-pnorm((p[2L]-p[3L])/p[4L]))*
              exp(-(y-p[3L])^2/2/((p[4L])^2+x^2))
      return(-sum(log(prop1+prop2)))
    } # end if
  } # end function minfunc
  ### function RmamErr to estimate parameters' Std.Err,
  ### added in 2013.07.25
  ### ***********************************************
    RmamErr<-function(pars, ED, sED, addsigma) {
      npars<-length(pars)
      nED<-length(ED)
      spars<-vector(length=npars)
      ### tol<-.Machine$double.eps^0.3
      message<-0
      ### finite-difference approximation hessian matrix
      res<-.Fortran("aperr",as.double(ED),as.double(sED),as.integer(nED),
            as.double(pars),as.double(addsigma),spars=as.double(spars),
            as.integer(npars),as.double(tol),message=as.integer(message),
            package="RadialPlotter")
      if(res$message!=0) {
        return(NULL)
      } else {
        return(res$spars)
      } # end if    
    } # end function RmamErr
    ### ***********************************************
    ### set boundaries for parameters
    if(ncomp==-1L) {
      lower<-c(1e-4, min(y), 1e-3)
      upper<-c(0.99, max(y), 5.0)
    } else if(ncomp==-2L) {
      lower<-c(1e-4, min(y), min(y), 1e-3)
      upper<-c(0.99, max(y), max(y), 5.0)
    } # end if
    ### find out three cluster using K-Means  
    ### then sorting them by ascending order
    kclus<-suppressWarnings(try(kmeans(x=y, centers=3L, iter.max=20L, nstart=100L), silent=TRUE))
    ### don't forget always doing a error checking
    if(class(kclus)=="try-error") stop("Error: k-means fails in parameter initializing!")
    kclus<-sort(kclus$centers)
    ### set gamma, mu and sdy values
    gama<-kclus[1L]
    mu<-c(kclus[2L], mean(kclus), mean(y))
    sdy<-sd(y)
    ### start the optimization by various
    ### initial values calling nlminb
    cmaxlik<-maxlik<-1e30
    errorflag<-1
    bexist<-0
    if(ncomp==-1L) {
      ### for three parameters (MAM3)
      for(i in 1:3) {
        for(j in 1:5) {
          ini<-c(0.01*(5.0)^(i-1),gama,sdy*0.4*j)
          res<-suppressWarnings(try(nlminb(start=ini,objective=minfunc,scale=1,
               control=list(iter.max=maxiter),lower=lower,upper=upper),silent=TRUE))
          if(class(res)!="try-error" && res$convergence==0) { 
            aperr<-RmamErr(res$par,drop(EDdata[,1L]),drop(EDdata[,2L]),addsigma=addsigma)
            if(!is.null(aperr) &&                 # std.error can be approximated
               res$objective<maxlik &&            # decreasing in minus maxlik
               all(abs(res$par-lower)>=1e-4) &&   # within low boundaries
               all(abs(res$par-upper)>=1e-4) ) {  # within up boundaries
              pars<-res$par
              error<-aperr
              maxlik<-res$objective
              errorflag<-0
            } # end if
            if(!is.null(aperr) &&                 # std.error can be approximated
               res$objective<cmaxlik &&           # decreasing in minus maxlik
               errorflag==1) {                    # might near to the boundaries
              cpars<-res$par
              cerror<-aperr
              cmaxlik<-res$objective
              bexist<-1
            } # end if
          } # end if
        } # end for
      } # end for
    } else if(ncomp==-2L) {
      ### for four parameters (MAM4)
      for(i in 1:3) {
        for(j in 1:3) {
          for(k in 1:3) {
            ini<-c(0.01*(5.0)^(i-1), gama, mu[k],sdy*0.5*j)
            res<-suppressWarnings(try(nlminb(start=ini,objective=minfunc,scale=1,
                 control=list(iter.max=maxiter),lower=lower,upper=upper),silent=TRUE))
            if(class(res)!="try-error" && res$convergence==0) { 
               aperr<-RmamErr(res$par,drop(EDdata[,1L]),drop(EDdata[,2L]),addsigma=addsigma)
               if(!is.null(aperr) &&                 # std.error can be approximated
                  res$objective<maxlik &&            # decreasing in minus maxlik
                  all(abs(res$par-lower)>=1e-4) &&   # within low boundaries
                  all(abs(res$par-upper)>=1e-4) &&   # within up boundaries
                  res$par[2L]<=res$par[3L] ) {       # gamma must not larger than mu
                 pars<-res$par
                 error<-aperr
                 maxlik<-res$objective
                 errorflag<-0
               } # end if
              if(!is.null(aperr) &&                  # std.error can be approximated
                  res$objective<cmaxlik &&           # decreasing in minus maxlik
                  res$par[2L]<=res$par[3L] &&        # gamma must not larger than mu
                  errorflag==1) {                    # might near to the boundaries
                 cpars<-res$par
                 cerror<-aperr
                 cmaxlik<-res$objective
                 bexist<-1
              } # end if
            } # end if
          } # end for
        } # end for
      } # end for
    } # end if
    ### now output the results
    if(errorflag==0) {
      ### reset parameters before output
      if(ncomp==-1L) {
        pars[2L]<-exp(pars[2L])
        error[2L]<-pars[2L]*error[2L]
      } else if(ncomp==-2L) {
        pars[2:3]<-exp(pars[2:3])
        error[2:3]<-pars[2:3]*error[2:3]
      } # end if
      return(list("errorflag"=0,"maxlik"=-maxlik,
                  "pars"=cbind(pars,error),
                  "BIC"=ifelse(ncomp==-1L,
                   2*maxlik+3*log(length(y)),
                   2*maxlik+4*log(length(y)))))
    } else if(bexist==1){
      ### reset parameters before output
      if(ncomp==-1L) {
        cpars[2L]<-exp(cpars[2L])
        cerror[2L]<-cpars[2L]*cerror[2L]
      } else if(ncomp==-2L) {
        cpars[2:3]<-exp(cpars[2:3])
        cerror[2:3]<-cpars[2:3]*cerror[2:3]
      } # end if
      return(list("errorflag"=1,"maxlik"=-cmaxlik,
                  "pars"=cbind(cpars,cerror),
                  "BIC"=ifelse(ncomp==-1L,
                   2*cmaxlik+3*log(length(y)),
                   2*cmaxlik+4*log(length(y)))))
    } else {
      return(NULL)
    } # end if
  } # end function Rmam
  ### ******************************************************
  ### now calculate the commom age model based equivalent dose
  z.i<-log(ed)
  se.i<-sqrt( (error/ed)^2 + addsigma^2 )
  commonED<-sum(z.i/se.i^2)/sum(1/se.i^2)
  scommonED <- 1/sqrt(sum(1/se.i^2))
  ### start do optimization
  errorflag<-0
  maxlik<-BIC<-0
  loopcomp<-0
  BILI<-NULL
  ### for CAM and FMM analysis
  if(ncomp%in%(0:maxcomp))  {
    if (ncomp==0L)  {
      goodcomp<-0
      BILI<-matrix(0, nrow=maxcomp-1, ncol=2L)
      ### call Fortran subroutine 'FineComp' to pick out the appropriate 
      ### number of component automatically
      Results<-.Fortran("FineComp",as.double(ed),as.double(error),as.integer(n),
                goodcomp=as.integer(goodcomp),as.double(addsigma),as.integer(maxiter),       
                as.double(eps),as.integer(maxcomp),BILI=as.double(BILI),pacakge="RadialPlotter")
      BILI<-Results$BILI
      dim(BILI)<-c(maxcomp-1, 2L)
      colnames(BILI)<-c("BIC","Maxlik")
      rownames(BILI)<-paste(rep("k=", maxcomp-1), c(2:maxcomp), sep ="") 
      ncomp<-Results$goodcomp
      loopcomp<-1
    } # end if
    spars<-pars<-matrix(0,nrow=2L, ncol=ncomp)
    ### call Fortran subroutine 'FineED1' to calculate parameters for 
    ### central age model or finite mixture age model
    Results<-.Fortran("FineED1",as.double(ed),as.double(error),as.integer(n),as.integer(ncomp),
              as.double(addsigma),pars=as.double(pars),spars=as.double(spars),maxlik=as.double(maxlik),
              BIC=as.double(BIC),as.integer(maxiter),as.double(eps),as.double(tol),
              errorflag=as.integer(errorflag),pacakge="RadialPlotter")  
    ### store errorflag, BIC, maxlik, ParsAndErrors
    errorflag<-Results$errorflag
    BIC<-Results$BIC
    maxlik<-Results$maxlik
    ParsAndErrors<-cbind(matrix(Results$pars, byrow=TRUE, ncol=2L),
                         matrix(Results$spars, byrow=TRUE, ncol=2L))
    if(ncomp==1L)  { 
      ### set matrix ParsAndErrors for CAM analysis
      ParsAndErrors<-matrix(ParsAndErrors[,c(1L,3L,2L,4L)])
      rownames(ParsAndErrors)<-c("Sigma", "Std.Sigma", "ED", "Std.ED")
      colnames(ParsAndErrors)<-"CAM.Value"
      ### plot or not (CAM)
      if(plot==TRUE)  {
        RadialPlot(Data=EDdata, Pars=ParsAndErrors[3L], zscale=zscale, 
                   samplename=samplename, kratio=kratio, pcolor=pcolor, psize=psize)
      } # end if
    }  else {  
      ### set matrix ParsAndErrors for FMM analysis
      ParsAndErrors<-ParsAndErrors[,c(1L,3L,2L,4L)][order(ParsAndErrors[,2L]), , drop=FALSE]
      colnames(ParsAndErrors)<-c("P","Std.P","ED","Std.ED") 
      rownames(ParsAndErrors) <- paste(rep("comp", ncomp), c(1:ncomp), sep = "") 
      ### plot or not (FMM)
      if(plot==TRUE)  {
        RadialPlot(Data=EDdata, Pars=ParsAndErrors[,3L], zscale=zscale, 
                   samplename=samplename, kratio=kratio, pcolor=pcolor, psize=psize)
      } # end if
    } # end if
  }  else  { 
    ### set default algorithm
    algorithm<-algorithm[1L]
    ### for MAM analysis 
    if(algorithm=="lbfgsb") {
      ### for algorithm 'lbfgsb'
      pars<-spars<-vector(length=2-ncomp)
      ### call Fortran subroutine 'MAM' to fit the minimum age models
      Results<-.Fortran("MAM",as.double(ed),as.double(error),as.integer(n),pars=as.double(pars),
                spars=as.double(spars), maxlik=as.double(maxlik), BIC=as.double(BIC),as.integer(2-ncomp),
                as.double(addsigma), as.integer(maxiter),as.double(tol),errorflag=as.integer(errorflag),
                pacakge="RadialPlotter")
      ### store errorflag, BIC, maxlik, ParsAndErrors
      errorflag<-Results$errorflag
      BIC<-Results$BIC
      maxlik<-Results$maxlik
      ParsAndErrors<-cbind(Results$pars, Results$spars)
    } else if(algorithm=="port") {
      ### for algorithm 'port'
      Results<-Rmam(EDdata, ncomp=ncomp, addsigma=addsigma, maxiter=maxiter)
      ### checking error of Rmam and store errorflag, BIC, maxlik, ParsAndErrors
      if(!is.null(Results)) {
        errorflag<-Results$errorflag
        BIC<-Results$BIC
        maxlik<-Results$maxlik
        ParsAndErrors<-Results$pars
      } else {
        ParsAndErrors<-matrix(0, nrow=2-ncomp, ncol=2L) 
        errorflag<-1
        BIC<-maxlik<-0
      } # end if
    } # end if
    ### set matrix ParsAndErrors for MAM
    colnames(ParsAndErrors)<-c("Pars", "Std.Pars")
    if(ncomp==-1L)  {
      rownames(ParsAndErrors)<-c("P", "gamma", "sigma")
    } else if(ncomp==-2L) {
      rownames(ParsAndErrors)<-c("P", "gamma", "mu", "sigma")
    } # end if
    ### plot or not (MAM)
    if(plot==TRUE)  {
      if(any(ParsAndErrors[,2L]==0)==FALSE)  {
        RadialPlot(Data=EDdata, Pars=ParsAndErrors[2L,1L], zscale=zscale, 
                   samplename=samplename, kratio=kratio, pcolor=pcolor, psize=psize)
      } else {
        ### alternatively, plot commonED out
        RadialPlot(Data=EDdata, Pars=exp(commonED), zscale=zscale, 
                   samplename=samplename, kratio=kratio, pcolor=pcolor, psize=psize)
      } # end if
    } # end if
  } # end if
  ### at this point the age model analysis have been terminated
  ### organize the returned list
  out<-list("errorflag"=errorflag,
            "commonED"=exp(commonED)*c(1,scommonED),
            "ncomp"=ncomp,       "maxcomp"=maxcomp, "loopcomp"=loopcomp,
            "pars"=ParsAndErrors,"BIC"=BIC,         "maxlik"=maxlik,"BILI"=BILI)
  ### set out to be of class 'RadialPlotter'
  class(out)<-"RadialPlotter"
  ### set out invisible
  invisible(out)
} # end function RadialPlotter.default
###
### set print method for object 'RadialPlotter', 2013.07.25
print.RadialPlotter<-function(x,...) {
  ### output the result to the terminal screen
  cat("\n")
  cat("==================Results of RadialPlotter==================","\n\n")
  ### error message
  cat("Error message:",x$errorflag,"\n\n")
  ### common ED
  cat("Common equivalent dose:",round(x$commonED[1L],5L),"+-",round(x$commonED[2L],5L),"\n\n")
  if(x$ncomp==1L)  {
    ### overdispersion in CAM
    cat("Overdispersion and Equivalent Dose are:","\n\n")
  }  else if(x$ncomp%in%c(0,2:x$maxcomp))  {
    ### FMM 
    if(x$loopcomp==1L)  {
      ### best number of components
      cat("The best component number is estimated at:",x$ncomp,"\n\n")
    } # end if
    ### parameters for CAM, FMM
    cat("Parameters and standard errors are:","\n")
  }  else  {
    ### MAM
    if(any(x$pars[,2L]==0)==FALSE)  {
      if(x$errorflag==1) {
        ### boundary checking
        cat("Warning: bad estimation, at least one parameter is near to the boundary!","\n\n")
      } # end if
      ### parameters for MAM
      cat("Parameters and standard errors are:","\n")
    }  else  {
      ### fails in MAM
      cat("Error  : parameters and standard errors can not be estimated!","\n")
      cat("Warning: commonED is drwan instead of MAMED!","\n\n")
    } # end if
  } # end if
  ### BIC, maxlik for CAM, FMM, MAM
  if(!x$ncomp %in% (-2:-1) || any(x$pars[,2L]==0)==FALSE)  {
    cat("----------------------------------","\n")
    print(round(as.data.frame(x$pars),5L))
    cat("----------------------------------","\n\n")
    cat("                      BIC value:",round(x$BIC,5L),"\n\n")
    cat("Maximum logged likelihood value:",round(x$maxlik,5L),"\n\n")
    ### looped BIC and maxlik values for FMM
    if(x$loopcomp==1L)  {
      cat("BIC values and Maximum logged likelihood values are:","\n")
      cat("------------------------","\n")
      print(round(as.data.frame(x$BILI),5L))
      cat("------------------------","\n\n")
    } # end if
  } # end if
  cat("===========================The End===========================","\n\n")
} # end function print.RadialPlotter
############################################ END FUNCTION RadialPlotter #######################################

################################################## FUNCTION calED #############################################
### *****************************************************************************************
### Function calED is used to calculate equivalent dose value.
###
###    Author: Peng Jun, 2013.06.22, revised in 2013.07.26, revised in 2013.08.01, revised in 2013.09.18
###
### Reference: Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived 
###            from single aliquot regenerative dose measurements. Ancient TL 25, pp. 15-24.
###            Duller, G., 2007. Analyst. pp. 27-28.
### *****************************************************************************************
calED<-
function(Curvedata,Ltx,model=c("line","exp","line+exp"),
         origin=FALSE,nstart=100,upb=1,ErrorMethod=c("mc","sp"),
         MCiter=1000,plot=TRUE,samplename=NULL) {
  UseMethod("calED")
}
### default method for function calED
calED.default<-
function(Curvedata,Ltx,model=c("line","exp","line+exp"),
         origin=FALSE,nstart=100,upb=1,ErrorMethod=c("mc","sp"),
         MCiter=1000,plot=TRUE,samplename=NULL)  {
  ### stop if not
  if(!is.data.frame(Curvedata))     stop("Error: Curvedata must be of type data.frame!")
  if(ncol(Curvedata)!=3L)           stop("Error: Curvedata must have three columns!")
  if(!is.numeric(
     unlist(unclass(Curvedata))))   stop("Error: elements in Curvedata must be all of type numeric!")
  if(any(Curvedata<0.0))            stop("Error: all elements in Curvedata must be not smaller than zero!")
  if(!is.numeric(Ltx))              stop("Error: Ltx must be a numeric vector!")
  if(length(Ltx)!=2L)               stop("Error: Ltx must contains two elements!")
  if(Ltx[1L]>max(Curvedata[,2L]))   stop("Error: Ltx must not exceed maximum Lx/Tx in Curvedata!")
  if(Ltx[1L]<=0.0)                  stop("Error: Ltx must larger than zero!")
  if(Ltx[2L]<=0.0)                  stop("Error: std.error of Ltx must larger than 0!")
  if(Ltx[2L]>max(Curvedata[,2L]))   stop("Error: std.error of Ltx must not exceed maximum Lx/Tx in Curvedata!")
  if(!is.character(model))          stop("Error: model must be 'line', 'exp' or 'line+exp'!")
  if(length(model)>=4L)             stop("Error: model must contain no more than three elements!")
  if(length(model)==1L) {
    if(!model %in% 
       c("line","exp","line+exp"))  stop("Error: model must be one of 'line', 'exp', 'line+exp'!")
  } else if(length(model)>=2L) {
    if(!all(model %in% 
       c("line","exp","line+exp"))) stop("Error: incorrect model!")
  } # end if
  if(!is.logical(origin))           stop("Error: origin must be of type logical!")
  if(length(origin)!=1L)            stop("Error: origin must be an one-element vector!")
  if(!is.numeric(nstart))           stop("Error: nstart must be of type numeric!")
  if(length(nstart)!=1L)            stop("Error: nstart must be an one-element vector!")
  if(abs(nstart-round(nstart))
     >=.Machine$double.eps^0.5)     stop("Error: nstart must be an integer!")
  if(nstart<=1L)                    stop("Error: nstart must larger than 1!")
  if(nstart>10000L)                 stop("Error: nstart is too large!")          
  if(!is.numeric(upb))              stop("Error: upb must be of type numeric!")
  if(length(upb)!=1L)               stop("Error: upb must be an one-element vector!")
  if(upb<=0)                        stop("Error: upb must larger than 0!")
  if(upb>1e3)                       stop("Error: upb is too large!")
  if(!is.character(ErrorMethod))    stop("Error: ErrorMethod must be of type character!")
  if(length(ErrorMethod)>=3L)       stop("Error: ErrorMethod must contain no more than two elements!")
  if(length(ErrorMethod)==1L) {
    if(!ErrorMethod %in%
       c("sp","mc"))                stop("Error: ErrorMethod must be either 'sp' or 'mc'!")
  } else if(length(ErrorMethod==2L)) {
    if(!all(ErrorMethod %in% 
       c("sp","mc")))               stop("Error: incorrect ErrorMethod!")
  } # end if
  if(!is.numeric(MCiter))           stop("Error: MCiter must be of type numeric!")
  if(length(MCiter)!=1L)            stop("Error: MCiter must be an one-element vector!")
  if(abs(MCiter-round(MCiter))
     >=.Machine$double.eps^0.5)     stop("Error: MCiter must be an integer!")
  if(MCiter<100L)                   stop("Error: MCiter is too small!")
  if(MCiter>50000L)                 stop("Error: MCiter is too large!")
  if(!is.logical(plot))             stop("Error: plot must be either TRUE or FALSE!")
  if(!is.character(samplename) &&
     !is.null(samplename))          stop("Error: samplename must be NULL or of type character!")
  ### set default model (linear)
  model<-model[1L]
  ### set default ErrorMethod (monte carlo)
  ErrorMethod<-ErrorMethod[1L]
  ### ************************************************
  ### check if data is enough for model fitting
  if(origin==TRUE) {
    if(model=="line" && length(levels(factor(Curvedata[,1L])))<1L) {
      stop("Error: fitting a linear model (origin) needs at least one paired independent observations")
    } # end if
  } else {
    if(model=="line" && length(levels(factor(Curvedata[,1L])))<2L) {
      stop("Error: fitting a linear model (non-origin) needs at least two paired independent observations")
    } # end if
  } # end if
  if(model=="exp" && length(levels(factor(Curvedata[,1L])))<3L) {
     stop("Error: fitting a exponential model needs at least three paired independent observations")
  } # end if
  if(model=='line+exp' && length(levels(factor(Curvedata[,1L])))<4L) {
    stop("fitting a linear+exponential model needs at least four paired independent observations")
  } # end if
  ### specify parameters that will be used 
  ### in Fortran subroutine 'calED.f90'
  Dose<-drop(Curvedata[,1L])                             # dose
  ndose<-nrow(Curvedata)                                 # ndose
  ltx<-cbind(drop(Curvedata[,2L]), drop(Curvedata[,3L])) # ltx
  inltx<-Ltx                                             # inltx
  outDose<-vector(length=2L)                             # outDose
  npars<-if(origin==TRUE) {
    if(model=="line") {
    1L}  else if(model=="exp") {
    2L} else if(model=="line+exp") {
    3L} # end if                           
    } else {
    if(model=="line") {
    2L}  else if(model=="exp") {
    3L} else if(model=="line+exp") {
    4L} # end if      
    } # end if                             # npars
  pars<-vector(length=npars)               # pars
  predtval<-vector(length=ndose)           # predtval
  parserrors<-vector(length=npars)         # parserrors
  value<- -99.0                            # value
  mcED<-vector(length=MCiter)              # mcED
  method<-ifelse(ErrorMethod=="sp",1L,2L)  # method
  motoiter<-MCiter                         # motoiter
  errorflag<-vector(length=2L)             # errorflag
  ### ************************************************
  ### call Fortran subroutine calED
  fFortran<-if(origin==FALSE) "calED" else "calED2"
  res<-.Fortran(fFortran,as.double(Dose),as.double(ltx),as.integer(ndose),as.double(inltx),
                outDose=as.double(outDose),pars=as.double(pars),as.integer(npars),predtval=as.double(predtval),
                as.double(upb),as.integer(nstart),parserrors=as.double(parserrors),value=as.double(value),
                mcED=as.double(mcED),as.integer(method),as.integer(motoiter),
                errorflag=as.integer(errorflag),package="RadialPlotter")
  ### error checking 
  if(res$errorflag[1L]!=123) {
    stop('Error: fail in Dose-Response curve fitting, fitting model might be inappropriate, or upb need to be modified!')
  } # end if
  ### set LMpars for output
  LMpars<-cbind(res$pars, res$parserrors)
  colnames(LMpars)<-c("Pars", "Std.Pars")
  rowname<-c("a", "b", "c", "d")
  rownames(LMpars)<-rowname[1:npars]   
  if(res$errorflag[2L]!=0)  {
    LMpars[,2L]<-NA
  } # end if
  ### set fit.value for output
  fit.value<-cbind(Curvedata[,c(1L,2L)], res$predtval)
  colnames(fit.value)<-c("ReDose", "Lx/Tx", "Fit.Lx/Tx")
  rownames(fit.value)<-paste("ReDose.", 1:ndose, sep="")
  ### set Dose for output
  ED<-c("ED"=res$outDose[1L], "Std.ED"=res$outDose[2L])
  ### reset mcED
  if(ErrorMethod=="sp") {
    res$mcED<-NULL
  } # end if
  ### prepare results for output
  output<-list("mcED"=res$mcED,      "LMpars"=LMpars,
               "residual"=res$value, "fit.value"=fit.value, "ED"=ED)
  ### plot or not
  if(plot==TRUE) {
    ### set pars
    par(bg="grey95",
        mgp=c(2,1,0),
        mar=c(3,3,2,1)+0.1)
    Xlim<-max(Curvedata[,1L], ED[1L])
    Ylim<-max(Curvedata[,2L], inltx[1L])
    plot(NA, NA, main=samplename, xlab="ReDose (Gy)",ylab="Lx/Tx",las=0,
         xlim=c(0,Xlim*1.05), ylim=c(0,Ylim*1.05), xaxs="i", yaxs="i", lab=c(7,7,9))
    ### add a filled density curve if ErrorMethod='mc'
    if(ErrorMethod=="mc") {
      dmcED<-density(res$mcED)
      dxy<-data.frame(unclass(dmcED)[c(1L,2L)])
      dxy[,2L]<-(dxy[,2L]-min(dxy[,2L]))/(max(dxy[,2L])-min(dxy[,2L]))*Ltx[1L]*0.9
      polygon(dxy, col="grey")
      rug(res$mcED)
    } # end if
    ### add ReDose as points
    points(Curvedata[,c(1L,2L)], pch=21, cex=3, bg="white")
    ### add error bars to Lx/Tx
    if(any(ltx[,2L]>=1e-3)) {
      arr<-suppressWarnings(try(arrows(Dose[ltx[,2L]>=1e-3], ltx[ltx[,2L]>=1e-3,1L]-ltx[ltx[,2L]>=1e-3,2L]/2, Dose[ltx[,2L]>=1e-3],
                                       ltx[ltx[,2L]>=1e-3,1L]+ltx[ltx[,2L]>=1e-3,2L]/2, code=3, lwd=2.5, angle=90, length=0.05, col="black"), 
                                       silent=TRUE))
    } # end if
    ### points calculate Equivalent Dose .VS. Ltx
    points(ED[1L], inltx[1L], pch=23, cex=3, bg="grey")
    ### add a fitting curve to the plot
    x<-NULL
    if(origin==FALSE) {
      if(npars==2L) {
        curve(LMpars[1L,1L]*x+LMpars[2L,1L], type="l",
              add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } else if(npars==3L) {
        curve(LMpars[1L,1L]*(1-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L], type="l",
              add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } else if(npars==4L) {
        curve(LMpars[1L,1L]*(1-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L]*x+LMpars[4L,1L],
              type="l", add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } # end if
    } else {
      if(npars==1L) {
        curve(LMpars[1L,1L]*x, type="l", add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } else if(npars==2L) {
        curve(LMpars[1L,1L]*(1-exp(-LMpars[2L,1L]*x)), type="l",
              add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } else if(npars==3L) {
        curve(LMpars[1L,1L]*(1-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L]*x,
              type="l", add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } # end if
    } # end if
    ### add a dash line to the plot
    lines(c(0,ED[1L],ED[1L]), c(inltx[1L],inltx[1L],0), lty="dashed", lwd=2)
    ### add error bars if standard errors for Doses are available
    if(ED[2L]>=1e-3) {
      arr<-suppressWarnings(try(arrows(ED[1L]-ED[2L]/2, inltx[1L], ED[1L]+ED[2L]/2,
                                inltx[1L], code=3, lwd=2.5, angle=90, length=0.05, col="black"), 
                                silent=TRUE))
    } # end if
    ### add a legend
    Curvetype<-if(model=="line") {
                 "Linear" } else if(model=="exp") {
                 "Exponential" } else if(model=="line+exp") {
                 "Exponential plus Linear"} # end if
    legend("topleft", legend=c(paste("Curve type: ", Curvetype,sep=""),
           paste("ED=", round(ED[1L],2L), " +- ", round(ED[2L],3L), " (Gy)", sep="")),
           yjust=2, ncol=1, cex=1.05, bty="o")
  grid()
  box(lwd=2)
  on.exit(par(bg="transparent",
              mgp=c(3,1,0),
              mar=c(5,4,4,2)+0.1))
  } # end if
  ### output the results    
  invisible(output)
} # end function calED
####################################### END FUNCTION calED ##################################################

######################################### FUNCTION decomp ###################################################
### ********************************************************************************************************
### Function decomp is used to decompose OSL signal
### curve, either type 'CW' or 'LM' can be analyzed.
###
###      Author: Peng Jun, 2013.06.05, revised in 2013.07.25, revised in 2013.09.18.
###
###  References: Bluszcz, A., 1996. Exponential function fitting to TL growth data and similar
###              applications. Geochronometria 13, 135–141.
###
###              Bluszcz, A., Adamiec, G., 2006. Application of differential evolution to fitting
###              OSL decay curves. Radiation Measurements 41, 886-891.
###
###              Jain, M., Murray, A.S., Boetter-Jensen, L., 2003. Characterisation of blue-light stimulated 
###              luminescence components in different quartz samples: implications for dose measurement. 
###              Radiation Measurements, 37 (4-5), pp. 441-449.
### ********************************************************************************************************
decomp<-
function(Sigdata,ncomp=3,typ=c("cw","lm"),
         LEDpower=60,LEDwavelength=470,
         plot=TRUE,xlog=TRUE,xylim=NULL,lwd=1,
         samplename=NULL,outfile=NULL,control.args=list()) {
  UseMethod("decomp")
} 
### default method for decomp
decomp.default<-
function(Sigdata,ncomp=3,typ=c("cw","lm"),
         LEDpower=60,LEDwavelength=470,
         plot=TRUE,xlog=TRUE,xylim=NULL,lwd=1,
         samplename=NULL,outfile=NULL,control.args=list()) {
  ### stop if not
  if(!is.data.frame(Sigdata))       stop("Error: Sigdata must be of type data.frame!")
  if(ncol(Sigdata)!=2L)             stop("Error: Sigdata must has two columns!")
  if(!is.numeric(
     unlist(unclass(Sigdata))))     stop("Error: elements in Sigdata must be all of type numeric!")
  if(any(Sigdata[,1L]<0))           stop("Error: the first column of Sigdata [time] must not contain values that below zero!")
  if(!is.numeric(ncomp))            stop("Error: ncomp must be of type numeric!")
  if(length(ncomp)!=1L)             stop("Error: ncomp must be an one-element vector!")
  if(!ncomp %in% (1:7) )            stop("Error: ncomp must be a integer ranges from 1 to 7!")
  if(!is.character(typ))            stop("Error: typ must be of type character!")
  if(length(typ)==1L) {
    if(!typ %in% c("cw","lm"))      stop("Error: typ must be either 'cw' or 'lm'!")
    if(typ=="cw" && 
       which.max(Sigdata[,2L])>3L)  stop("Error: incorrect typ!")
    if(typ=="lm" &&
       which.max(Sigdata[,2L])<=3L) stop("Error: incorrect typ!")
  } else if(length(typ)==2L) {
    if(!all(typ %in% c("cw","lm"))) stop('Error: incorrect typ!')
    if(typ[1L]=="cw" && 
       which.max(Sigdata[,2L])>3L)  stop("Error: incorrect typ!")
    if(typ[1L]=="lm" &&
       which.max(Sigdata[,2L])<=3L) stop("Error: incorrect typ!")
  } else {                          stop("Error: typ must contains one or two elements!")
  } # end if
  if(!is.numeric(LEDpower))         stop("Error: LEDpower must be of type numeric!")
  if(length(LEDpower)!=1L)          stop("Error: LEDpower must be an one-element vector!")
  if(LEDpower<=0.0)                 stop("Error: LEDpower must be larger than zero!")
  if(!is.numeric(LEDwavelength))    stop("Error: LEDwavelength must be of type numeric!")
  if(length(LEDwavelength)!=1L)     stop("Error: LEDwavelength must be an one-element vector!")
  if(LEDwavelength<=0.0)            stop("Error: LEDwavelength must be larger than zero!")
  if(!is.logical(plot))             stop("Error: plot must be either TRUE or FALSE!")
  if(length(plot)!=1L)              stop("Error: plot must be an one-element vector!")
  if(!is.logical(xlog))             stop("Error: xlog must be either TRUE or FALSE!")
  if(length(xlog)!=1L)              stop("Error: xlog must be an one-element vector!")
  if(!is.null(xylim) &&
     !is.numeric(xylim))            stop("Error: xylim must be NULL or of type numeric!")
  if(!is.null(xylim)) {
    if(length(xylim)!=2L)           stop("Error: xylim must be of length two!")
    if(any(xylim<=0.0))             stop("Error: xylim must be larger than zero!")
    if(xylim[1L]>max(Sigdata[,1L])) stop("Error: xylim[1] must not exceed the maximum stimulative time value!")
    if(xylim[2L]>max(Sigdata[,2L])) stop("Error: xylim[2] must not exceed the maximum signal value!")
  } # end if
  if(!is.numeric(lwd))              stop("Error: lwd must be of type numeric!")
  if(length(lwd)!=1L)               stop("Error: lwd must be an one-element vector!")
  if(!is.null(samplename) &&
     !is.character(samplename))     stop("Error: samplename must be NULL or type of character!")
  if(!is.null(samplename)) {
    if(length(samplename)!=1L)      stop("Error: sample name must be an one-element vector!")
  } # end if
  if(!is.null(outfile) &&
     !is.character(outfile))        stop("Error: outfile must be NULL or type of character!")
  if(!is.null(outfile)) {
    if(length(outfile)!=1L)         stop("Error: outfile must be an one-element vector!")
  } # end if
  if(class(control.args)!="list")   stop("Error: control.args must be a list!")
  if(length(control.args)>=6L)      stop("Error: number of parameters in control.args must smaller than 6!")
  if(!all(names(control.args) %in% 
     list("factor","f","cr","maxiter","tol"))) stop("Error: incorrect parameter in control.args!")
  ### public parameters for subroutine  
  ### 'decomp' and subroutine 'fitlm'
  tim<-drop(Sigdata[,1L])
  sig<-drop(Sigdata[,2L])
  ntim<-length(tim)
  pars<-Stdpars<-vector(length=2*ncomp)
  value<--99.0
  predtval<-vector(length=ntim)
  ### set default typ (CW)
  typ<-typ[1L]
  ### do type dependent fitting
  if(typ=="cw") {
    ### for OSL signal of type 'CW'
    ### default arguments for differential evolution
    ### ****************************************************
    args<-list(factor=10L,f=0.5,cr=0.99,maxiter=1000L,tol=0.1)   
    args[names(control.args)]<-control.args
    factor<-args$factor
    f<-args$f
    cr<-args$cr
    maxiter<-args$maxiter
    tol<-args$tol
    if(!is.numeric(factor))       stop("Error: control.args: factor must be of type numeric!")
    if(length(factor)!=1L)        stop("Error: control.args: factor must be an one-element vector!")
    if(abs(factor-round(factor))
       >=.Machine$double.eps^0.5) stop("Error: control.args: factor must be an integer!")
    if(factor<3L)                 stop("Error: control.args: factor is too small!") 
    if(factor>20L)                stop("Error: control.args: factor is too large!")             
    if(!is.numeric(f))            stop("Error: control.args: f must be of type numeric!")
    if(length(f)!=1L)             stop("Error: control.args: f must be an one-element vector!")
    if(f<=0.0 || f>1.2)           stop("Error: control.args: incorrect f!")
    if(!is.numeric(cr))           stop("Error: control.args: cr must be of type numeric!")
    if(length(cr)!=1L)            stop("Error: control.args: cr must be an one-element vector!")
    if(cr<=0.0 || cr>1.0)         stop("Error: control.args: incorrect cr!")
    if(!is.numeric(maxiter))      stop("Error: control.args: maxiter must be of type numeric!")
    if(length(maxiter)!=1L)       stop("Error: control.args: maxiter must be an one-element vector!")
    if(abs(maxiter-round(maxiter))
       >=.Machine$double.eps^0.5) stop("Error: control.args: maxiter must be an integer!")
    if(maxiter<10L)               stop("Error: control.args: maxiter is too small!")
    if(maxiter>10000L)            stop("Error: control.args: maxiter is too large!") 
    if(!is.numeric(tol))          stop("Error: control.args: tol must be of type numeric!")
    if(length(tol)!=1L)           stop("Error: control.args: tol must be an one-element vector!")
    if(tol<0.0)                   stop("Error: control.args: incorrect tol!")
    ### ****************************************************
    errorflag<-vector(length=3L)
    ### call subroutine decomp
    res<-.Fortran("decomp",as.integer(ncomp),as.double(tim),as.double(sig),
          as.integer(ntim),pars=as.double(pars),Stdpars=as.double(Stdpars),
          value=as.double(value),predtval=as.double(predtval),as.integer(factor),
          as.double(f),as.double(cr),as.integer(maxiter),as.double(tol),
          errorflag=as.integer(errorflag),package="RadialPlotter")
    ### Error checking
    if(res$errorflag[1L]==1 && 
       res$errorflag[3L]==-1)  {
       stop(paste("Error: signal can not be decomposed to",ncomp,"components!"))
    } # end if
  } else if(typ=="lm") {
    ### for OSL signal of type "LM"
    errorflag<-0
    res<-.Fortran("fitlm",as.integer(ncomp),as.double(tim),as.double(sig),as.integer(ntim),
          pars=as.double(pars),Stdpars=as.double(Stdpars),value=as.double(value),
          predtval=as.double(predtval),errorflag=as.integer(errorflag),
          package="RadialPlotter")
    ### error checking
    if(res$errorflag==1) {
      stop(paste("Error: signal can not be decomposed to",ncomp,"components!"))
    } # end if
  } # end if
  ### decide Photoionisation cross-section (cm2)
  h<-6.62606957e-34
  ny<-299792458/(LEDwavelength/10^9)
  E<-h*ny
  LEDpower<-LEDpower/1000
  ### reshpae pars for output
  pars<-cbind(res$pars[1:ncomp], res$Stdpars[1:ncomp], 
              res$pars[-(1:ncomp)], res$Stdpars[-(1:ncomp)],
              res$pars[-(1:ncomp)]/LEDpower*E)
  pars<-pars[order(pars[,3L], decreasing=TRUE), , drop=FALSE]
  colnames(pars)<-c("Ithn", "Std.Ithn", "Lamda", "Std.Lamda", "Pcs")
  rownames(pars)<-paste("Comp.", 1:ncomp, sep="")
  ### reset pars and errorflag
  if(typ=="cw") {
    if(res$errorflag[3L]==-1) {
      pars[,c(2L,4L)]<-NA
      errorflag<-1
    } else {
      errorflag<-0
    } # end if
  } else if(typ=="lm") {
    if(res$errorflag==-1) {
      pars[,c(2L,4L)]<-NA
      errorflag<-1
    } else {
      errorflag<-0
    } # end if
  } # end if
  ### Calculate signal values for each component
  CompSig<-apply(cbind(pars[,1L], pars[,3L]), MARGIN=1,
           function(x) if(typ=="cw") x[1L]*exp(-x[2L]*tim) else 
           x[1L]*x[2L]*(tim/max(tim))*exp(-x[2L]*tim^2/2/max(tim)))
  SigProp<-colSums(CompSig)/sum(sig)
  CompSig<-round(cbind(res$predtval, CompSig), 5L)
  colnames(CompSig)<-c("Fit.Signal", paste("Comp.", 1:ncomp, sep=""))
  ### plot or not
  if(plot==TRUE) {
    ### set pars
    par(bg="grey95",
        mgp=c(2,1,0),
        mar=c(3,3,2,1)+0.1)
    ### reset tim, sig, CompSig for ploting accroding to xylim
    if(!is.null(xylim)) {
      indexts<-tim<=xylim[1L]
      tim<-tim[indexts]
      sig<-sig[indexts]
      CompSig1<-CompSig[indexts, ,drop=FALSE]
    } else {
      CompSig1<-CompSig
    } # end if
    ### add a scatter plot (Time .VS. Signal)
    plot(tim, sig, main=samplename, log=ifelse(xlog==TRUE,"x",""), las=0, lab=c(7,7,9),
         ylim=c(-max(sig)*0.01, ifelse(is.null(xylim), max(sig), min(max(sig), xylim[2L]))*1.01),
         xlab="Stimulated Time (s)", ylab="OSL Counts", xaxs="r", yaxs="i", 
         type="p", pch=21, cex=1, bg="White", col="black") 
    ### set colors
    colors<-c("blue", "red", "green", "yellow", "purple", "orange", "brown")
    ### lines Time .VS. Fitted values
    lines(tim, CompSig1[,1L], lwd=lwd, col="black", lty="solid")
    ### lines Time .VS. Component signal (1 to ncomp)
    for(i in 1:(ncomp-1)) {
        lines(tim, CompSig1[,i+1], lwd=lwd, col=colors[i], lty="solid")
      } # end for
    ### lines Time .VS. Back.ground signal
    lines(tim, CompSig1[,ncomp+1], lwd=lwd, col=colors[7L], lty="solid")
    ### add a legend
    legend(ifelse(typ=="cw", "topright", ifelse(xlog==TRUE, "topleft", ifelse(which.max(sig)>length(sig)/2, "topleft","topright"))), 
           legend=c("Fitted.Curve", paste("Comp.", 1:ncomp," (",round(SigProp*100,2L),"%)", sep="")),
           col=c("black", colors[c(1:(ncomp-1),7)]), pch=c(21,rep(NA,ncomp)), lty="solid",
           yjust=2, ncol=1, cex=1, bty="o", lwd=lwd, pt.bg="White")
  } # end if
  grid(equilogs=FALSE)
  box(lwd=2)
  ### reset pars
  on.exit(par(bg="transparent",
              mgp=c(3,1,0),
              mar=c(5,4,4,2)+0.1))
  ### results for output
  out<-list("Comp.Signal"=CompSig, "pars"=pars,
            "value"=res$value,     "errorflag"=errorflag)
  ### wirte fitted signal values to a file or not
  if(!is.null(outfile)) {
    write.csv(CompSig, file=paste(outfile, ".csv"))
  } # end if
  ### out is invisible
  invisible(out)
} # end function decomp
######################################## END FUNCTION decomp ####################################

########################################## FUNCTION sgcED #######################################
### ***********************************************************************************************
### Function sgcED is used to analyze equivalent doses using SGC method.
###
###     Author: Peng Jun, 2013.06.23, revised in 2013.07.26, revised in 2013.08.02, revised in 2013.09.18
###
### References: Roberts,H.M. and Duller,G.A.T., 2004. Standardised growth curves for optical dating 
###             of sediment using multiple-grain aliquots. Radiation Measurements 38, pp. 241-252.
###
###             Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from 
###             single aliquot regenerative dose measurements. Ancient TL 25, pp. 15-24.
### ***********************************************************************************************
sgcED<-
function(Curvedata,Ltx,model=c("line","exp","line+exp"),
         origin=FALSE,nstart=100,upb=1,ErrorMethod=c("mc","sp"),
         MCiter=1000,plot=TRUE,samplename=NULL,outfile=NULL)  {
   UseMethod("sgcED")
}
### default method for function sgc.ED
sgcED.default<-
function(Curvedata,Ltx,model=c("line","exp","line+exp"),
         origin=FALSE,nstart=100,upb=1,ErrorMethod=c("mc","sp"),
         MCiter=1000,plot=TRUE,samplename=NULL,outfile=NULL)  {
  ### stop if not
  if(!is.data.frame(Curvedata))          stop("Error: Curvedata must be of type data.frame!")
  if(ncol(Curvedata)!=3L)                stop("Error: Curvedata must have three columns!")
  if(!is.numeric(
     unlist(unclass(Curvedata))))        stop("Error: elements in Curvedata must be all of type numeric!")
  if(any(Curvedata<0.0))                 stop("Error: all elements in Curvedata must be not smaller than zero!")
  if(!is.data.frame(Ltx))                stop("Error: Ltx must be a data.frame!")
  if(ncol(Ltx)!=2L)                      stop("Error: Ltx must contains two columns!")
  if(!is.numeric(
     unlist(unclass(Ltx))))              stop("Error: elements in Ltx must be all of type numeric!")
  if(any(Ltx[,1L]>max(Curvedata[,2L])))  stop("Error: Ltx must not exceed maximum Lx/Tx in Curvedata!")
  if(any(Ltx[,1L]<=0.0))                 stop("Error: Ltx must larger than zero!")
  if(any(Ltx[,2L]<=0.0))                 stop("Error: std.error of Ltx must larger than zero!")
  if(any(Ltx[,2L]>max(Curvedata[,2L])))  stop("Error: std.error of Ltx must not exceed maximum Lx/Tx in Curvedata!")
  if(!is.character(model))               stop("Error: model must be of type character!")
  if(length(model)>=4L)                  stop("Error: model must contain no more than three elements!")
  if(length(model)==1L) {
    if(!model %in% 
       c("line","exp","line+exp"))       stop("Error: model must be one of 'line', 'exp', 'line+exp'!")
  } else if(length(model)>=2L) {
    if(!all(model %in% 
       c("line","exp","line+exp")))      stop("Error: incorrect model!")
  } # end if
  if(!is.logical(origin))                stop("Error: origin must be of type logical!")
  if(length(origin)!=1L)                 stop("Error: origin must be an one-element vector!")
  if(!is.numeric(nstart))                stop("Error: nstart must be of type numeric!")
  if(length(nstart)!=1L)                 stop("Error: nstart must be an one-element vector!")
  if(abs(nstart-round(nstart))
     >=.Machine$double.eps^0.5)          stop("Error: nstart must be an integer!")
  if(nstart<=1L)                         stop("Error: nstart must larger than 1!")
  if(nstart>10000L)                      stop("Error: nstart is too large!")          
  if(!is.numeric(upb))                   stop("Error: upb must be of type numeric!")
  if(length(upb)!=1L)                    stop("Error: upb must be an one-element vector!")
  if(upb<=0)                             stop("Error: upb must larger than 0!")
  if(upb>1e3)                            stop("Error: upb is too large!")
  if(!is.character(ErrorMethod))         stop("Error: ErrorMethod must be of type character!")
  if(length(ErrorMethod)>=3L)            stop("Error: ErrorMethod must contain no more than two elements!")
  if(length(ErrorMethod)==1L) {
    if(!ErrorMethod %in%
       c("sp","mc"))                     stop("Error: ErrorMethod must be either 'sp' or 'mc'!")
  } else if(length(ErrorMethod==2L)) {
    if(!all(ErrorMethod %in% 
       c("sp","mc")))                    stop("Error: incorrect ErrorMethod!")
  } # end if
  if(!is.numeric(MCiter))                stop("Error: MCiter must be of type numeric!")
  if(length(MCiter)!=1L)                 stop("Error: MCiter must be an one-element vector!")
  if(abs(MCiter-round(MCiter))
     >=.Machine$double.eps^0.5)          stop("Error: MCiter must be an integer!")
  if(MCiter<100L)                        stop("Error: MCiter is too small!")
  if(MCiter>10000L)                      stop("Error: MCiter is too large!")
  if(!is.logical(plot))                  stop("Error: plot must be either TRUE or FALSE!")
  if(!is.character(samplename) &&
     !is.null(samplename))               stop("Error: samplename must be NULL or of type character!")
  if(!is.null(samplename)) {
    if(length(samplename)!=1L)           stop("Error: samplename must be an one-element vector!")
  } # end if
  if(!is.null(outfile) &&
     !is.character(outfile))             stop("Error: outfile must be NULL or of type character!")
  if(!is.null(outfile)) {
    if(length(outfile)!=1L)              stop("Error: outfile must be an one-element vector!")
  } # end if
  ### set default model (linear)
  model<-model[1L]
  ### set default ErrorMethod ('mc')
  ErrorMethod<-ErrorMethod[1L]
  ### check if data is enough for model fitting
  if(origin==TRUE) {
    if(model=="line" && length(levels(factor(Curvedata[,1L])))<1L) {
      stop("Error: fitting a linear model (origin) needs at least one paired independent observations")
    } # end if
  } else {
    if(model=="line" && length(levels(factor(Curvedata[,1L])))<2L) {
      stop("Error: fitting a linear model (non-origin) needs at least two paired independent observations")
    } # end if
  } # end if
  if(model=="exp" && length(levels(factor(Curvedata[,1L])))<3L) {
     stop("Error: fitting a exponential model needs at least three paired independent observations")
  } # end if
  if(model=="line+exp" && length(levels(factor(Curvedata[,1L])))<4L) {
    stop("fitting a linear+exponential model needs at least four paired independent observations")
  } # end if
  ### parameters for subroutine 'sgcED.f90'
  Dose<-drop(Curvedata[,1L])
  ltx<-cbind(drop(Curvedata[,2L]), drop(Curvedata[,3L]))
  ndose<-length(Dose)
  inltx<-cbind(drop(Ltx[,1L]), drop(Ltx[,2L]))
  ninltx<-nrow(Ltx)
  outDose<-matrix(0,nrow=ninltx, ncol=2L)
  npars<-if(origin==TRUE) {
    if(model=="line") {
    1L } else if(model=="exp") {
    2L } else if(model=="line+exp") {
    3L } # end if
    } else {
    if(model=="line") {
    2L } else if(model=="exp") {
    3L } else if(model=="line+exp") {
    4L } # end if
  } # end if
  pars<-parserrors<-vector(length=npars)
  predtval<-vector(length=ndose)
  value<- -99.0
  method<-ifelse(ErrorMethod=="sp",1L,2L)
  errorflag<-vector(length=2L)
  ### calculate equivalent doses
  fFortran<-if(origin==FALSE) "sgcED" else "sgcED2"
  res<-.Fortran(fFortran,as.double(Dose),as.double(ltx),as.integer(ndose),as.integer(ninltx),
                as.double(inltx),outDose=as.double(outDose),pars=as.double(pars),as.integer(npars),predtval=as.double(predtval),
                as.double(upb),as.integer(nstart),parserrors=as.double(parserrors),value=as.double(value),
                as.integer(method),as.integer(MCiter),errorflag=as.integer(errorflag),package="RadialPlotter")
  ### error checking
  if(res$errorflag[1L]!=123) {
    stop("Error: fail in Dose-Response curve fitting, fitting model might be inappropriate, or upb need to be modified!")
  } # end if
  ### set LMpars
  LMpars<-cbind(res$pars, res$parserrors)
  colnames(LMpars)<-c("Pars", 'Std.Pars')
  rowname<-c("a", "b", "c", "d")
  rownames(LMpars)<-rowname[1:npars]    
  if(res$errorflag[2L]!=0) {
    LMpars[,2L]<-NA
  } # end if
  ### set fit.value for output
  fit.value<-cbind(drop(Curvedata[,1L]), drop(Curvedata[,2L]), res$predtval)
  colnames(fit.value)<-c("ReDose", "Lx/Tx", "Fit.Lx/Tx")
  rownames(fit.value)<-paste("ReDose.", 1:ndose, sep="")
  ### set ED values
  ED<-matrix(res$outDose,ncol=2L)
  rownames(ED)<-paste("NO.", 1:ninltx, sep="")
  colnames(ED)<-c("ED", "Std.ED")
  ### prepare results for output
  output<-list("LMpars"=LMpars,       "residual"=res$value,
               "fit.value"=fit.value, "ED"=ED)
  ### plot or not
  if(plot==TRUE) {
    ### set pars
    par(bg="grey95",
        mgp=c(2,1,0),
        mar=c(3,3,2,1)+0.1)
    ### set limitations for x and y axies
    Xlim<-max(Curvedata[,1L],ED[,1L])
    Ylim<-max(Curvedata[,2L],Ltx[,1L])
    ### plot a dose-response curve
    plot(NA, NA, main=samplename, las=0,
         xlab="ReDose (Gy)", ylab="Lx/Tx",
         xlim=c(0,Xlim*1.05), ylim=c(0,Ylim*1.05),
         xaxs="i", yaxs="i", lab=c(7,7,9) )
    ### add ReDose as points
    points(Curvedata[,c(1L,2L)], pch=21, cex=3, bg="white")
    ### points calculate Equivalent Dose .VS. Ltx
    points(ED[,1L], Ltx[,1L], pch=23, cex=3, bg='grey')
    ### add error bars to Lx/Tx
    if(any(Curvedata[,3L]>=1e-3)) {
      arr<-suppressWarnings(try(arrows(drop(Curvedata[Curvedata[,3L]>=1e-3,1L]), 
                                       drop(Curvedata[Curvedata[,3L]>=1e-3,2L])-
                                       drop(Curvedata[Curvedata[,3L]>=1e-3,3L])/2, 
                                       drop(Curvedata[Curvedata[,3L]>=1e-3,1L]),
                                       drop(Curvedata[Curvedata[,3L]>=1e-3,2L])+
                                       drop(Curvedata[Curvedata[,3L]>=1e-3,3L])/2, 
                                       code=3, lwd=2.5, angle=90, length=0.05, col="black"), silent=TRUE))
    } # end if
    ### add error bars to calculated ED
    if(any(ED[,2L]>=1e-3)) {
      arr<-suppressWarnings(try(arrows(ED[ED[,2L]>=1e-3,1L]-ED[ED[,2L]>=1e-3,2L]/2,Ltx[ED[,2L]>=1e-3,1L],
                                       ED[ED[,2L]>=1e-3,1L]+ED[ED[,2L]>=1e-3,2L]/3, Ltx[ED[,2L]>=1e-3,1L], 
                                       code=3, lwd=2.5, angle=90, length=0.05, col="black"), silent=TRUE))
    } # end if
    ### reset model
    model<-model[1L]
    ### add a fitting curve
    x<-NULL
    if(origin==TRUE)  {
      if(model=="line") {
        curve(LMpars[1L,1L]*x, type="l",
              add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } else if(model=="exp") {
        curve(LMpars[1L,1L]*(1-exp(-LMpars[2L,1L]*x)), type="l",
              add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } else if(model=="line+exp")  {
        curve(LMpars[1L,1L]*(1-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L]*x, type="l",
              add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } # end if
    } else {
      if(model=="line") {
        curve(LMpars[1L,1L]*x+LMpars[2L,1L], type="l",
              add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } else if(model=="exp") {
        curve(LMpars[1L,1L]*(1-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L], type="l",
              add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } else if(model=="line+exp")  {
        curve(LMpars[1L,1L]*(1-exp(-LMpars[2L,1L]*x))+LMpars[3L,1L]*x+LMpars[4L,1L], type="l",
              add=TRUE, lw=2, from=0, to=Xlim*1.05)
      } # end if
    } # end if
    ### add dash lines 
    for(i in 1: ninltx) {
      lines(c(0,ED[i,1L],ED[i,1L]),c(Ltx[i,1L],Ltx[i,1L],0),lty="dashed",lwd=1)
    } # end for
  grid()
  box(lwd=2)
  ### reset pars
  on.exit(par(bg="transparent",
              mgp=c(3,1,0),
              mar=c(5,4,4,2)+0.1))
  } # end if
  ### write equivalent doses out or not
  if(!is.null(outfile)) {
    write.csv(ED,file=paste(outfile,".csv"))
  } # end if
  return(output)
} # end function sgcED   
######################################## END FUNCTION sgcED ###############################################
