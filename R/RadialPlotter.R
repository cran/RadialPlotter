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
###     Author: Peng Jun, 2013.04.20, revised in 2013.07.26.
###
### References: Galbraith, R.F., 1988. Graphical Display of Estimates Having Differing Standard 
###             Errors. Technometrics, 30 (3), pp. 271-281.
###
###             Galbraith, R.F., Roberts, R.G., 2012. Statistical aspects of equivalent dose and error calculation and  
###             display in OSL dating: An overview and some recommendations. Quaternary Geochronology, 11, pp. 1-27.
### ******************************************************************************************************************
RadialPlotter<-
function(EDdata,ncomp=0,addsigma=0,maxiter=500,
         maxcomp=9,algorithm=c('lbfgsb','port'),
         eps=.Machine$double.eps^0.5,plot=TRUE,
         zscale=NULL,samplename=NULL)  {
  UseMethod('RadialPlotter')
}
### set default method for function RadialPlotter
RadialPlotter.default<-
function(EDdata,ncomp=0,addsigma=0,maxiter=500,
         maxcomp=9,algorithm=c('lbfgsb','port'),
         eps=.Machine$double.eps^0.5,plot=TRUE,
         zscale=NULL,samplename=NULL)  {
  ### stop if not 
  if(!is.data.frame(EDdata))    stop('Error: EDdata must be of type data.frame!')
  if(ncol(EDdata)!=2L)          stop('Error: EDdata must contain two columns!')
  if(!is.numeric(
     unlist(unclass(EDdata))))  stop('Error: elements in EDdata must be all of type numeric!')
  if(any(EDdata[,1]<=0))        stop('Error: equivalent dose must be larger than 0!')
  if(!is.numeric(maxcomp))      stop('Error: maxcomp must be of type numeric!')
  if(length(maxcomp)!=1L)       stop('Error: maxcomp must be an one-element vector!')
  if(abs(maxcomp-round(maxcomp))
     >=.Machine$double.eps^0.5) stop('Error: maxcomp must be a integer!')                                
  if(maxcomp<0)                 stop('Error: maxcomp must not be smaller than 0!')
  if(maxcomp>nrow(EDdata))      stop('Error: maxcomp must not exceed the length of EDdata!')
  if(!is.numeric(ncomp))        stop('Error: ncomp must be of type numeric!')
  if(length(ncomp)!=1L)         stop('Error: ncomp must be an one-element vector!')
  if(!ncomp %in% (-2:maxcomp))  stop('Error: ncomp must be an integer ranging from -2 to maxcomp!')
  if(!is.numeric(addsigma))     stop('Error: addsigma must be of type numeric!')
  if(length(addsigma)!=1L)      stop('Error: addsigma must be an one-element vector!')
  if(addsigma<0)                stop('Error: addsigma must not be smaller than 0!')
  if(!is.numeric(maxiter))      stop('Error: maxiter must be of type numeric!')
  if(length(maxiter)!=1L)       stop('Error: maxiter must be an one-element vector!')
  if(abs(maxiter-round(maxiter))
     >=.Machine$double.eps^0.5) stop('Error: maxiter must be an integer!')
  if(maxiter<10L)               stop('Error: maxiter is too small!')
  if(maxiter>10000L)            stop('Error: maxiter is too large!')
  if(!is.character(algorithm))  stop('Error: algorithm must be of type character!')
  if(length(algorithm)==1L) {
    if(!algorithm %in% c('lbfgsb','port'))      stop('Error: algorithm must be either "lbfgsb" or "port"!')
  } else if(length(algorithm)==2L) {
    if(!all(algorithm %in% c('lbfgsb','port'))) stop('Error: only algorithm "lbfgsb" or "port" is available!')
  } else {                                      stop('Error: algorithm must contain no more than two elements!')
  } # end if
  if(!is.numeric(eps))          stop('Error: eps must be of type numeric!')
  if(length(eps)!=1L)           stop('Error: eps must be an one-element vector!')
  if(eps<0)                     stop('Error: eps must be not smaller than 0!')
  if(!is.logical(plot))         stop('Error: plot must be either TRUE or FALSE!')
  if(!is.null(zscale) &&
     !is.numeric(zscale))       stop('Error: zscale must either be NULL or a numeric vector!')
  if(!is.null(samplename) &&
     !is.character(samplename)) stop('Error: samplename must either be NULL or of type character!')
  if(ncomp==0L && maxcomp<=1L)  stop("Error: maxcomp must non't be smaller than 2 if ncomp is 0!")  
  ### set n, ED, sED
  n<-nrow(EDdata) 
  ed<-drop(EDdata[,1])
  error<-drop(EDdata[,2])
  ### function for radial plot drawing,
  ### the code is revised from Rex Galbraith 
  ### *************************************************
  RadialPlot<-function(Data,Pars,zscale,samplename)  {
    z.i<-log(Data[,1])
    se.i<-Data[,2]/Data[,1]
    ### if Pars=NULL, lines will not be plotted out
    if(!is.null(Pars))  {
      Pars<-log(Pars)
    } # end if
    plot.area_ratio<-4.5/6
    z.0<-sum(z.i/se.i^2)/sum(1/se.i^2)
    x.i<-1/se.i
    y.i<-(z.i-z.0)/se.i
    h<-plot.area_ratio*(max(x.i)-min(x.i))/(max(y.i)-min(y.i))
    r.0<-max(sqrt(x.i^2+y.i^2*h^2))
    if(is.null(zscale))  {
      mkest<-round(exp(pretty(z.i)),0L)
      if(sd(z.i)/mean(z.i)<0.05)  {
        mkest<-round(exp(pretty(c(min(z.i)*0.98,max(z.i)*1.02))),0L)
      } # end if
    } else  {
      mkest<-zscale
    } # end if
    mkr<-range(mkest)
    circ.x<-(0:200)/200
    circ.x<-mkr[1]+circ.x*(mkr[2]-mkr[1])
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
    par(mfrow=c(1,1),oma=c(2,0,0,0),mar=c(5,4,4,4.5),xpd=TRUE,las=1,cex=1)
    plot(NA,NA,xlim=c(0,max(x.i)),ylim=c(-yaxis.max,yaxis.max),xlab='',ylab='',
       xaxt='n',yaxt='n',main=samplename,bty='n',typ='n',xaxs='i',yaxs='i')      
    lenpar<-length(Pars) 
    maxx<-which.max(circ.x)
    ### if Pars!=NULL, starting drawing lines out
    if(!is.null(Pars))  {
      for(i in 1:lenpar)  {
        if(Pars[i]-z.0>=0)  {
          Ci<-approx(circ.y[maxx:201],circ.x[maxx:201],ties="ordered",
          xout=(Pars[i]-z.0)*max(x.i),rule=2,method="constant",f=0.5)$y
        }  else  {
          Ci<-approx(circ.y[1:(maxx-1)],circ.x[1:(maxx-1)],ties="ordered",
          xout=(Pars[i]-z.0)*max(x.i),rule=2,method="constant",f=0.5)$y
        } # end if
       lines(c(0,Ci),c(0,(Pars[i]-z.0)*Ci),lty=1,col='black',lwd=1.5)
      } # end for
    } # end if
    lines(circ.x,circ.y,lwd=2)
    segments(xmk1,ymk1,xmk2,ymk2)
    text(xmk3,ymk3,round(mkest,2L),cex=1)
    text(max(x.i)*1.15,0,'De (Gy)',cex=1,srt=90)
    par(mfrow=c(1,1),oma=c(2,0,0,0.5),mar=c(5,4,4,4.5),xpd=TRUE,las=1,cex=1,new=TRUE)
    plot(NA,NA,xlim=c(0,max(x.i)),ylim=c(-yaxis.max, yaxis.max),yaxt='n',xaxt='n',
    xaxs='i',yaxs='i',ylab='Standardised Estimate',xlab='',typ='n',bty='n',)
    points(x.i,y.i,pch=1,col='black',cex=1.5)
    axis(side=1,line=4,cex.axis=1,lwd=2)
    mtext(side=1,line=6,'Precision',cex=1)
    reticks.labels<-round(1/axTicks(side=1)*100,1L)
    reticks.values<-axTicks(side=1)
    axis(side=1,at=reticks.values[-1], lwd=2,labels=reticks.labels[-1],line=4,cex.axis=1,tck=0.02,padj=-4)
    mtext(side=1,line=1.5,'Relative Error (%)',cex=1)
    axis(side=2,at=c(-2,-1,0,1,2),lwd=2,labels=c('-2','','0','','2'),cex.axis=1)
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
    x<-sqrt((EDdata[,2]/EDdata[,1])^2+addsigma^2)
    ### change Equivalent Dose to log-scale
    y<-log(drop(EDdata[,1]))
    ### function minfunc to be minimized
    ### ******************************
    minfunc<-function(p)  {
    ### here ncomp is a gloable variable
    if(ncomp==-1L) {
      ### for MAM3
      u0<-(p[2]/(p[3])^2+y/x^2)/(1/(p[3])^2+1/x^2)
      sigma0<-1/(sqrt(1/(p[3])^2+1/x^2)) 
      prop1<-p[1]/sqrt(2*pi)/x*exp(-(y-p[2])^2/2/x^2)
      prop2<-(1-p[1])/sqrt(2*pi*((p[3])^2+x^2))* 
             (1-pnorm((p[2]-u0)/sigma0))/(1-pnorm((p[2]-p[2])/p[3]))*
              exp(-(y-p[2])^2/2/((p[3])^2+x^2))
      return(-sum(log(prop1+prop2)))
    } else if(ncomp==-2L) {
      ### for MAM4
      u0<-(p[3]/(p[4])^2+y/x^2)/(1/(p[4])^2+1/x^2)
      sigma0<-1/(sqrt(1/(p[4])^2+1/x^2)) 
      prop1<-p[1]/sqrt(2*pi)/x*exp(-(y-p[2])^2/2/x^2)
      prop2<-(1-p[1])/sqrt(2*pi*((p[4])^2+x^2))* 
             (1-pnorm((p[2]-u0)/sigma0))/(1-pnorm((p[2]-p[3])/p[4]))*
              exp(-(y-p[3])^2/2/((p[4])^2+x^2))
      return(-sum(log(prop1+prop2)))
    } # end if
  } # end function minfunc
  ### function RmamErr to estimate parameters' Std.Err,
  ### added in 2013.07.25
  ### ***********************************************
    RmamErr<-function(pars,ED,sED,addsigma) {
      npars<-length(pars)
      nED<-length(ED)
      spars<-vector(length=npars)
      ### tol<-.Machine$double.eps^0.3
      message<-0
      ### finite-difference approximation hessian matrix
      res<-.Fortran('aperr',as.double(ED),as.double(sED),as.integer(nED),
            as.double(pars),as.double(addsigma),spars=as.double(spars),
            as.integer(npars),as.double(tol),message=as.integer(message),
            package='RadialPlotter')
      if(res$message!=0) {
        return(NULL)
      } else {
        return(res$spars)
      } # end if    
    } # end function RmamErr
    ### ***********************************************
    ### set boundaries for parameters
    if(ncomp==-1L) {
      lower<-c(1e-4,min(y),1e-3)
      upper<-c(0.99,max(y),5.0)
    } else if(ncomp==-2L) {
      lower<-c(1e-4,min(y),min(y),1e-3)
      upper<-c(0.99,max(y),max(y),5.0)
    } # end if
    ### find out three cluster using K-Means  
    ### then sorting them by ascending order
    kclus<-suppressWarnings(try(kmeans(x=y,centers=3,iter.max=20,nstart=100),silent=TRUE))
    ### don't forget always doing a error checking
    if(class(kclus)=='try-error') stop('Error: k-means fails in parameter initializing!')
    kclus<-sort(kclus$centers)
    ### set gamma, mu and sdy values
    gama<-kclus[1]
    mu<-c(kclus[2],mean(kclus),mean(y))
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
          if(class(res)!='try-error' && res$convergence==0) { 
            aperr<-RmamErr(res$par,drop(EDdata[,1]),drop(EDdata[,2]),addsigma=addsigma)
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
            if(class(res)!='try-error' && res$convergence==0) { 
               aperr<-RmamErr(res$par,drop(EDdata[,1]),drop(EDdata[,2]),addsigma=addsigma)
               if(!is.null(aperr) &&                 # std.error can be approximated
                  res$objective<maxlik &&            # decreasing in minus maxlik
                  all(abs(res$par-lower)>=1e-4) &&   # within low boundaries
                  all(abs(res$par-upper)>=1e-4) &&   # within up boundaries
                  res$par[2]<=res$par[3] ) {         # gamma must not larger than mu
                 pars<-res$par
                 error<-aperr
                 maxlik<-res$objective
                 errorflag<-0
               } # end if
              if(!is.null(aperr) &&                  # std.error can be approximated
                  res$objective<cmaxlik &&           # decreasing in minus maxlik
                  res$par[2]<=res$par[3] &&          # gamma must not larger than mu
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
        pars[2]<-exp(pars[2])
        error[2]<-pars[2]*error[2]
      } else if(ncomp==-2L) {
        pars[2:3]<-exp(pars[2:3])
        error[2:3]<-pars[2:3]*error[2:3]
      } # end if
      return(list('errorflag'=0,'maxlik'=-maxlik,
                  'pars'=cbind(pars,error),
                  'BIC'=ifelse(ncomp==-1L,
                   2*maxlik+3*log(length(y)),
                   2*maxlik+4*log(length(y)))))
    } else if(bexist==1){
      ### reset parameters before output
      if(ncomp==-1L) {
        cpars[2]<-exp(cpars[2])
        cerror[2]<-cpars[2]*cerror[2]
      } else if(ncomp==-2L) {
        cpars[2:3]<-exp(cpars[2:3])
        cerror[2:3]<-cpars[2:3]*cerror[2:3]
      } # end if
      return(list('errorflag'=1,'maxlik'=-cmaxlik,
                  'pars'=cbind(cpars,cerror),
                  'BIC'=ifelse(ncomp==-1L,
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
      BILI<-matrix(0,nrow=maxcomp-1,ncol=2L)
      ### call Fortran subroutine 'FineComp' to pick out the appropriate 
      ### number of component automatically
      Results<-.Fortran('FineComp',as.double(ed),as.double(error),as.integer(n),
                goodcomp=as.integer(goodcomp),as.double(addsigma),as.integer(maxiter),       
                as.double(eps),as.integer(maxcomp),BILI=as.double(BILI),pacakge='RadialPlotter')
      BILI<-Results$BILI
      dim(BILI)<-c(maxcomp-1,2L)
      colnames(BILI)<-c('BIC','Maxlik')
      rownames(BILI)<-paste(rep("k=", maxcomp-1),c(2:maxcomp),sep ="") 
      ncomp<-Results$goodcomp
      loopcomp<-1
    } # end if
    spars<-pars<-matrix(0,nrow=2L,ncol=ncomp)
    ### call Fortran subroutine 'FineED1' to calculate parameters for 
    ### central age model or finite mixture age model
    Results<-.Fortran('FineED1',as.double(ed),as.double(error),as.integer(n),as.integer(ncomp),
              as.double(addsigma),pars=as.double(pars),spars=as.double(spars),maxlik=as.double(maxlik),
              BIC=as.double(BIC),as.integer(maxiter),as.double(eps),as.double(tol),
              errorflag=as.integer(errorflag),pacakge='RadialPlotter')  
    ### store errorflag, BIC, maxlik, ParsAndErrors
    errorflag<-Results$errorflag
    BIC<-Results$BIC
    maxlik<-Results$maxlik
    ParsAndErrors<-cbind(matrix(Results$pars,byrow=TRUE,ncol=2L),
                         matrix(Results$spars,byrow=TRUE,ncol=2L))
    if(ncomp==1L)  { 
      ### set matrix ParsAndErrors for CAM analysis
      ParsAndErrors<-matrix(ParsAndErrors[,c(1,3,2,4)])
      rownames(ParsAndErrors)<-c('Sigma','sSigma','ED','sED')
      colnames(ParsAndErrors)<-'CAM'
      ### plot or not (CAM)
      if(plot==TRUE)  {
        RadialPlot(EDdata,ParsAndErrors[3],zscale,samplename)
      } # end if
    }  else {  
      ### set matrix ParsAndErrors for FMM analysis
      ParsAndErrors<-ParsAndErrors[,c(1,3,2,4)][order(ParsAndErrors[,2]),]
      colnames(ParsAndErrors)<-c('P','sP','ED','sED') 
      rownames(ParsAndErrors) <- paste(rep("comp", ncomp), c(1:ncomp), sep = "") 
      ### plot or not (FMM)
      if(plot==TRUE)  {
        RadialPlot(EDdata,ParsAndErrors[,3],zscale,samplename)
      } # end if
    } # end if
  }  else  { 
    ### set default algorithm
    algorithm<-algorithm[1]
    ### for MAM analysis 
    if(algorithm=='lbfgsb') {
      ### for algorithm 'lbfgsb'
      pars<-spars<-vector(length=2-ncomp)
      ### call Fortran subroutine 'MAM' to fit the minimum age models
      Results<-.Fortran('MAM',as.double(ed),as.double(error),as.integer(n),pars=as.double(pars),
                spars=as.double(spars), maxlik=as.double(maxlik), BIC=as.double(BIC),as.integer(2-ncomp),
                as.double(addsigma), as.integer(maxiter),as.double(tol),errorflag=as.integer(errorflag),
                pacakge='RadialPlotter')
      ### store errorflag, BIC, maxlik, ParsAndErrors
      errorflag<-Results$errorflag
      BIC<-Results$BIC
      maxlik<-Results$maxlik
      ParsAndErrors<-cbind(Results$pars,Results$spars)
    } else if(algorithm=='port') {
      ### for algorithm 'port'
      Results<-Rmam(EDdata,ncomp=ncomp,addsigma=addsigma,maxiter=maxiter)
      ### checking error of Rmam and store errorflag, BIC, maxlik, ParsAndErrors
      if(!is.null(Results)) {
        errorflag<-Results$errorflag
        BIC<-Results$BIC
        maxlik<-Results$maxlik
        ParsAndErrors<-Results$pars
      } else {
        ParsAndErrors<-matrix(0,nrow=2-ncomp,ncol=2L) 
        errorflag<-1
        BIC<-maxlik<-0
      } # end if
    } # end if
    ### set matrix ParsAndErrors for MAM
    colnames(ParsAndErrors)<-c('pars','error')
    if(ncomp==-1L)  {
      rownames(ParsAndErrors)<-c('P','gama','sigma')
    } else if(ncomp==-2L) {
      rownames(ParsAndErrors)<-c('P','gama','mu','sigma')
    } # end if
    ### plot or not (MAM)
    if(plot==TRUE)  {
      if(any(ParsAndErrors[,2]==0)==FALSE)  {
        RadialPlot(EDdata,ParsAndErrors[2,1],zscale=zscale,samplename=samplename)
      } else {
        ### alternatively, plot commonED out
        RadialPlot(EDdata,Pars=exp(commonED),zscale=zscale,samplename=samplename)
      } # end if
    } # end if
  } # end if
  ### at this point the age model analysis have been terminated
  ### organize the returned list
  out<-list('errorflag'=errorflag,
            'commonED'=exp(commonED)*c(1,scommonED),
            'ncomp'=ncomp,'maxcomp'=maxcomp,'loopcomp'=loopcomp,
            'pars'=ParsAndErrors,'BIC'=BIC,'maxlik'=maxlik,'BILI'=BILI)
  ### set out to be of class 'RadialPlotter'
  class(out)<-'RadialPlotter'
  ### set out invisible
  invisible(out)
} # end function RadialPlotter.default
### set print method for object 'RadialPlotter', 2013.07.25
print.RadialPlotter<-function(x,...) {
  ### output the result to the terminal screen
  cat('\n')
  cat('==================Results for RadialPlotter==================','\n\n')
  ### error message
  cat('Error message:',x$errorflag,'\n\n')
  ### common ED
  cat('Common equivalent dose:',round(x$commonED[1],5L),'+-',round(x$commonED[2],5L),'\n\n')
  if(x$ncomp==1L)  {
    ### overdispersion in CAM
    cat('Overdispersion and Equivalent Dose are:','\n\n')
  }  else if(x$ncomp%in%c(0,2:x$maxcomp))  {
    ### FMM 
    if(x$loopcomp==1L)  {
      ### best number of components
      cat('The best component number is estimated at:',x$ncomp,'\n\n')
    } # end if
    ### parameters for CAM, FMM
    cat('Parameters and standard errors are:','\n')
  }  else  {
    ### MAM
    if(any(x$pars[,2]==0)==FALSE)  {
      if(x$errorflag==1) {
        ### boundary checking
        cat('Warning: bad estimation, at least one parameter is near to the boundary!','\n\n')
      } # end if
      ### parameters for MAM
      cat('Parameters and standard errors are:','\n')
    }  else  {
      ### fails in MAM
      cat('Error  : parameters and standard errors can not be estimated!','\n')
      cat('Warning: commonED is drwan instead of MAMED!','\n\n')
    } # end if
  } # end if
  ### BIC, maxlik for CAM, FMM, MAM
  if(!x$ncomp %in% (-2:-1) || any(x$pars[,2]==0)==FALSE)  {
    cat('----------------------------------','\n')
    print(round(as.data.frame(x$pars),5L))
    cat('----------------------------------','\n\n')
    cat('                      BIC value:',round(x$BIC,5L),'\n\n')
    cat('Maximum logged likelihood value:',round(x$maxlik,5L),'\n\n')
    ### looped BIC and maxlik values for FMM
    if(x$loopcomp==1L)  {
      cat('BICs and Maximum logged likelihood values are:','\n')
      cat('------------------------','\n')
      print(round(as.data.frame(x$BILI),5L))
      cat('------------------------','\n\n')
    } # end if
  } # end if
  cat('===========================The End===========================','\n\n')
} # end function print.RadialPlotter
############################################ END FUNCTION RadialPlotter #######################################
################################################## FUNCTION calED #############################################
### *****************************************************************************************
### Function calED is used to calculate equivalent dose value.
###
###    Author: Peng Jun, 2013.06.22, revised in 2013.07.26, revised in 2013.08.01
###
### Reference: Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived 
###            from single aliquot regenerative dose measurements. Ancient TL 25, pp. 15-24.
###            Duller, G., 2007. Analyst. pp. 27-28.
### *****************************************************************************************
calED<-
function(Curvedata,Ltx,model=c('line','exp','line+exp'),
         nstart=100,upb=5,ErrorMethod=c('mc','sp'),
         MCiter=1000,plot=TRUE,samplename=NULL) {
  UseMethod('calED')
}
### default method for function calED
calED.default<-
function(Curvedata,Ltx,model=c('line','exp','line+exp'),
         nstart=100,upb=5,ErrorMethod=c('mc','sp'),
         MCiter=1000,plot=TRUE,samplename=NULL)  {
  ### stop if not
  if(!is.data.frame(Curvedata))     stop('Error: Curvedata must be of type data.frame!')
  if(ncol(Curvedata)!=3L)           stop('Error: Curvedata must have three columns!')
  if(!is.numeric(
     unlist(unclass(Curvedata))))   stop('Error: elements in Curvedata must be all of type numeric!')
  if(!is.numeric(Ltx))              stop('Error: Ltx must be a numeric vector!')
  if(length(Ltx)!=2L)               stop('Error: Ltx must contains two elements!')
  if(any(Ltx<=0))                   stop('Error: Ltx must larger than 0!')
  if(any(Ltx>=max(Curvedata[,2])))  stop('Error: Ltx must not beyond maximum Lx/Tx in Curvedata!')
  if(!is.character(model))          stop('Error: model must be of type character!')
  if(length(model)>=4L)             stop('Error: model must contain no more than three elements!')
  if(length(model)==1L) {
    if(!model %in% 
       c('line','exp','line+exp'))  stop('Error: model must be one of "line", "exp", "line+exp"!')
  } else if(length(model)>=2L) {
    if(!all(model %in% 
       c('line','exp','line+exp'))) stop('Error: incorrect model!')
  } # end if
  if(!is.numeric(nstart))           stop('Error: nstart must be of type numeric!')
  if(length(nstart)!=1L)            stop('Error: nstart must be an one-element vector!')
  if(abs(nstart-round(nstart))
     >=.Machine$double.eps^0.5)     stop('Error: nstart must be an integer!')
  if(nstart<=1L)                    stop('Error: nstart must larger than 1!')
  if(nstart>10000L)                 stop('Error: nstart is too large!')          
  if(!is.numeric(upb))              stop('Error: upb must be of type numeric!')
  if(length(upb)!=1L)               stop('Error: upb must be an one-element vector!')
  if(upb<=0)                        stop('Error: upb must larger than 0!')
  if(upb>1e4)                       stop('Error: upb is too large!')
  if(!is.character(ErrorMethod))    stop('Error: ErrorMethod must be of type character!')
  if(length(ErrorMethod)>=3L)       stop('Error: ErrorMethod must contain no more than two elements!')
  if(length(ErrorMethod)==1L) {
    if(!ErrorMethod %in%
       c('sp','mc'))                stop('Error: ErrorMethod must be either "sp" or "mc"!')
  } else if(length(ErrorMethod==2L)) {
    if(!all(ErrorMethod %in% 
       c('sp','mc')))               stop('Error: incorrect ErrorMethod!')
  } # end if
  if(!is.numeric(MCiter))           stop('Error: MCiter must be of type numeric!')
  if(length(MCiter)!=1L)            stop('Error: MCiter must be an one-element vector!')
  if(abs(MCiter-round(MCiter))
     >=.Machine$double.eps^0.5)     stop('Error: MCiter must be an integer!')
  if(MCiter<100L)                   stop('Error: MCiter is too small!')
  if(MCiter>50000L)                 stop('Error: MCiter is too large!')
  if(!is.logical(plot))             stop('Error: plot must be either TRUE or FALSE!')
  if(!is.character(samplename) &&
     !is.null(samplename))          stop('Error: samplename must be NULL or of type character!')
  ### set default model (linear)
  model<-model[1]
  ### set default ErrorMethod (monte carlo)
  ErrorMethod<-ErrorMethod[1]
  ### check if data is enough for model fitting
  if(model=='line' && nrow(Curvedata)<2L) {
    stop('Error: fitting a linear model needs at least two paired observations')
  } # end if
  if(model=='exp' && nrow(Curvedata)<3L) {
     stop('Error: fitting a exponential model needs at least three paired observations')
  } # end if
  if(model=='line+exp' && nrow(Curvedata)<4L) {
    stop('fitting a linear+exponential model needs at least four paired observations')
  } # end if
  ### ************************************************
  ### specify parameters that will be used 
  ### in Fortran subroutine 'calED.f90'
  Dose<-drop(Curvedata[,1])                           # dose
  ndose<-nrow(Curvedata)                              # ndose
  ltx<-cbind(drop(Curvedata[,2]),drop(Curvedata[,3])) # ltx
  inltx<-Ltx                                          # inltx
  outDose<-vector(length=2L)                          # outDose
  npars<-if(model=='line') {
    2}  else if(model=='exp') {
    3} else if(model=='line+exp') {
    4} # end if                            # npars
  pars<-vector(length=npars)               # pars
  predtval<-vector(length=ndose)           # predtval
  parserrors<-vector(length=npars)         # parserrors
  value<- -99.0                            # value
  mcED<-vector(length=MCiter)              # mcED
  method<-ifelse(ErrorMethod=='sp',1,2)    # method
  motoiter<-MCiter                         # motoiter
  errorflag<-vector(length=2L)             # errorflag
  ### ************************************************
  ### call Fortran subroutine calED
  res<-.Fortran('calED',as.double(Dose),as.double(ltx),as.integer(ndose),as.double(inltx),
        outDose=as.double(outDose),pars=as.double(pars),as.integer(npars),predtval=as.double(predtval),
        as.double(upb),as.integer(nstart),parserrors=as.double(parserrors),value=as.double(value),
        mcED=as.double(mcED),as.integer(method),as.integer(motoiter),
        errorflag=as.integer(errorflag),package='RadialPlotter')
  ### error checking 
  if(res$errorflag[1]!=123) {
    stop('Error: fail in Dose-Response curve fitting, the fitting model might be inappropriate!')
  } # end if
  ### set LMpars for output
  LMpars<-cbind(res$pars,res$parserrors)
  colnames(LMpars)<-c('Estimate','Std.Error')
  rowname<-c('a','b','c','d')
  rownames(LMpars)<-rowname[1:npars]   
  if(res$errorflag[2]!=0) {
    LMpars[,2]<-NA
  } # end if 
  ### set fit.value for output
  fit.value<-cbind(Curvedata[,1:2],res$predtval)
  colnames(fit.value)<-c('ReDose','Lx/Tx','Fit.Lx/Tx')
  rownames(fit.value)<-paste('ReDose.',1:ndose,sep='')
  ### set Dose for output
  ED<-c('ED'=res$outDose[1],'Std.ED'=res$outDose[2])
  ### reset mcED
  if(ErrorMethod=='sp') {
    res$mcED<-NULL
  } # end if
  ### prepare results for output
  output<-list('mcED'=res$mcED,'LMpars'=LMpars,
               'residual'=res$value,'fit.value'=fit.value,'ED'=ED )
  ### plot or not
  if(plot==TRUE) {
    Xlim<-max(Curvedata[,1],ED[1])
    Ylim<-max(Curvedata[,2],inltx[1])
    plot(NA,NA,main=samplename,xlab='ReDose (Gy)',ylab='Lx/Tx',
         xlim=c(0,Xlim*1.05),ylim=c(0,Ylim*1.05),xaxs='i',yaxs='i',lab=c(7,7,9))
    ### add a filled density curve if ErrorMethod='mc'
    if(ErrorMethod=='mc') {
      dmcED<-density(res$mcED)
      dxy<-data.frame(unclass(dmcED)[1:2])
      dxy[,2]<-(dxy[,2]-min(dxy[,2]))/(max(dxy[,2])-min(dxy[,2]))*Ltx[1]*0.9
      polygon(dxy,col='grey')
    } # end if
    ### add ReDose as points
    points(Curvedata[,1:2],pch=1,cex=3)
    ### add error bars to Lx/Tx
    if(all(ltx[,2]>=1e-3)) {
      arr<-suppressWarnings(try(arrows(Dose,ltx[,1]-ltx[,2]/2,Dose,
           ltx[,1]+ltx[,2]/2,code=3,lwd=2.5,angle=90,length=0.05,col="black"),silent=TRUE))
    } # end if
    ### points calculate Equivalent Dose .VS. Ltx
    points(ED[1],inltx[1],pch=23,cex=3,bg='grey')
    ### add a fitting curve to the plot
    x<-NULL
    if(npars==2L) {
      curve(LMpars[1,1]*x+LMpars[2,1],type='l',
            add=TRUE,lw=2,from=0,to=Xlim*1.05)
    } else if(npars==3L) {
      curve(LMpars[1,1]*(1-exp(-LMpars[2,1]*x))+LMpars[3,1],type='l',
            add=TRUE,lw=2,from=0,to=Xlim*1.05)
    } else if(npars==4L) {
      curve(LMpars[1,1]*(1-exp(-LMpars[2,1]*x))+LMpars[3,1]*x+LMpars[4,1],type='l',
            add=TRUE,lw=2,from=0,to=Xlim*1.05)
    } # end if
    ### add a dash line to the plot
    lines(c(0,ED[1],ED[1]),c(inltx[1],inltx[1],0),lty='dashed',lwd=2)
    ### add error bars if standard errors for Doses are available
    if(ED[2]>=1e-3) {
      arr<-suppressWarnings(try(arrows(ED[1]-ED[2]/2,inltx[1],ED[1]+ED[2]/2,
           inltx[1],code=3,lwd=2.5,angle=90,length=0.05,col="black"),silent=TRUE))
    } # end if
    ### add a legend
    Curvetype<-if(model=='line') {
                 'Linear' } else if(model=='exp') {
                 'Exponential' } else if(model=='line+exp') {
                 'Exponential plus Linear'} # end if
    legend("topleft",legend=c(paste('Curve type: ',Curvetype,sep=''),
           paste('ED=',round(ED[1],3L),' +- ',round(ED[2],3L),' (Gy)',sep='')),
           yjust=2,ncol=1,cex=1.05,bty="o")
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
###      Author: Peng Jun, 2013.06.05, revised in 2013.07.25.
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
function(Sigdata,ncomp=3,typ=c('cw','lm'),
         LEDpower=60,LEDwavelength=470,
         plot=TRUE,lwd=1,xylim=NULL,samplename=NULL,
         outfile=NULL,control.args=list()) {
  UseMethod('decomp')
} 
### default method for decomp
decomp.default<-
function(Sigdata,ncomp=3,typ=c('cw','lm'),
         LEDpower=60,LEDwavelength=470,
         plot=TRUE,lwd=1,xylim=NULL,samplename=NULL,
         outfile=NULL,control.args=list()) {
  ### stop if not
  if(!is.data.frame(Sigdata))       stop('Error: Sigdata must be of type data.frame!')
  if(ncol(Sigdata)!=2L)             stop('Error: Sigdata must has two columns!')
  if(!is.numeric(
     unlist(unclass(Sigdata))))     stop('Error: elements in Sigdata must be all of type numeric!')
  if(any(Sigdata[,1]<0))            stop('Error: the first column of Sigdata (time) must not contain minus value!')
  if(!is.numeric(ncomp))            stop('Error: ncomp must be of type numeric!')
  if(length(ncomp)!=1L)             stop('Error: ncomp must be an one-element vector!')
  if(!ncomp %in% (1:7) )            stop('Error: ncomp must be a integer ranges from 1 to 7!')
  if(!is.character(typ))            stop('Error: typ must be of type character!')
  if(length(typ)==1L) {
    if(!typ %in% c('cw','lm'))      stop('Error: typ must be either "cw" or "lm"!')
  } else if(length(typ)==2L) {
    if(!all(typ %in% c('cw','lm'))) stop('Error: incorrect typ!')
  } else {                          stop('Error: typ must contains one or two elements!')
  } # end if
  if(!is.numeric(LEDpower))         stop('Error: LEDpower must be of type numeric!')
  if(length(LEDpower)!=1L)          stop('Error: LEDpower must be an one-element vector!')
  if(!is.numeric(LEDwavelength))    stop('Error: LEDwavelength must be of type numeric!')
  if(length(LEDwavelength)!=1L)     stop('Error: LEDwavelength must be an one-element vector!')
  if(!is.logical(plot))             stop('Error: plot must be either TRUE or FALSE!')
  if(!is.numeric(lwd))              stop('Error: lwd must be of type numeric!')
  if(length(lwd)!=1L)               stop('Error: lwd must be an one-element vector!')
  if(!is.null(xylim) && 
     !is.numeric(xylim))            stop('Error: xylim must be either NULL or a numeric vector!')
  if(!is.null(xylim)) {
    if(length(xylim)!=2L)           stop('Error: xylim must contain two elements!')
    if(any(xylim<=0))               stop('Error: xylim must larger than 0!')
    if(xylim[1]>
       max(drop(Sigdata[,1])))      stop('Error: upper boundary of X axis must not exceed maximum time of stimulation!')
  } # end if
  if(!is.null(samplename) &&
     !is.character(samplename))     stop('Error: samplename must be NULL or type of character!')
  if(!is.null(outfile) &&
     !is.character(outfile))        stop('Error: outfile must be NULL or type of character!')
  if(class(control.args)!='list')   stop('Error: control.args must be a list!')
  if(length(control.args)>=6L)      stop('Error: number of parameters in control.args must smaller than 6!')
  if(!all(names(control.args) %in% 
     list('factor','f','cr','maxiter','tol'))) stop('Error: incorrect parameter in control.args!')
  ### public parameters for subroutine  
  ### 'decomp' and subroutine 'fitlm'
  tim<-drop(Sigdata[,1])
  sig<-drop(Sigdata[,2])
  ntim<-length(tim)
  pars<-Stdpars<-vector(length=2*ncomp)
  value<--99.0
  predtval<-vector(length=ntim)
  ### set default typ (CW)
  typ<-typ[1]
  ### do type dependent fitting
  if(typ=='cw') {
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
    if(!is.numeric(factor))       stop('Error: control.args: factor must be of type numeric!')
    if(length(factor)!=1L)        stop('Error: control.args: factor must be an one-element vector!')
    if(abs(factor-round(factor))
       >=.Machine$double.eps^0.5) stop('Error: control.args: factor must be an integer!')
    if(factor<3L)                 stop('Error: control.args: factor is too small!') 
    if(factor>20L)                stop('Error: control.args: factor is too large!')             
    if(!is.numeric(f))            stop('Error: control.args: f must be of type numeric!')
    if(length(f)!=1L)             stop('Error: control.args: f must be an one-element vector!')
    if(f<=0.0 || f>1.2)           stop('Error: control.args: incorrect f!')
    if(!is.numeric(cr))           stop('Error: control.args: cr must be of type numeric!')
    if(length(cr)!=1L)            stop('Error: control.args: cr must be an one-element vector!')
    if(cr<=0.0 || cr>1.0)         stop('Error: control.args: incorrect cr!')
    if(!is.numeric(maxiter))      stop('Error: control.args: maxiter must be of type numeric!')
    if(length(maxiter)!=1L)       stop('Error: control.args: maxiter must be an one-element vector!')
    if(abs(maxiter-round(maxiter))
       >=.Machine$double.eps^0.5) stop('Error: control.args: maxiter must be an integer!')
    if(maxiter<10L)               stop('Error: control.args: maxiter is too small!')
    if(maxiter>10000L)            stop('Error: control.args: maxiter is too large!') 
    if(!is.numeric(tol))          stop('Error: control.args: tol must be of type numeric!')
    if(length(tol)!=1L)           stop('Error: control.args: tol must be an one-element vector!')
    if(tol<0.0)                   stop('Error: control.args: incorrect tol!')
    ### ****************************************************
    errorflag<-vector(length=3L)
    ### call subroutine decomp
    res<-.Fortran('decomp',as.integer(ncomp),as.double(tim),as.double(sig),
          as.integer(ntim),pars=as.double(pars),Stdpars=as.double(Stdpars),
          value=as.double(value),predtval=as.double(predtval),as.integer(factor),
          as.double(f),as.double(cr),as.integer(maxiter),as.double(tol),
          errorflag=as.integer(errorflag),package='RadialPlotter')
    ### Error checking
    if(res$errorflag[1]==1 && 
       res$errorflag[3]==-1)  {
       stop(paste('Error: signal can not be decomposed to',ncomp,'components!'))
    } # end if
  } else if(typ=='lm') {
    ### for OSL signal of type 'LM'
    errorflag<-0
    res<-.Fortran('fitlm',as.integer(ncomp),as.double(tim),as.double(sig),as.integer(ntim),
          pars=as.double(pars),Stdpars=as.double(Stdpars),value=as.double(value),
          predtval=as.double(predtval),errorflag=as.integer(errorflag),
          package='RadialPlotter')
    ### error checking
    if(res$errorflag==1) {
      stop(paste('Error: signal can not be decomposed to',ncomp,'components!'))
    } # end if
  } # end if
  ### decide Photoionisation cross-section (cm2)
  h<-6.62606957e-34
  ny<-299792458/(LEDwavelength/10^9)
  E<-h*ny
  LEDpower<-LEDpower/1000
  ### reshpae pars for output
  pars<-cbind(res$pars[1:ncomp],res$Stdpars[1:ncomp], 
              res$pars[-(1:ncomp)],res$Stdpars[-(1:ncomp)],
              res$pars[-(1:ncomp)]/LEDpower*E)
  pars<-pars[order(pars[,3],decreasing=TRUE),,drop=FALSE]
  colnames(pars)<-c('Ithn','Std.Ithn','Lamda','Std.Lamda','Pcs')
  rownames(pars)<-paste('Comp.',1:ncomp,sep='')
  ### reset pars and errorflag
  if(typ=='cw') {
    if(res$errorflag[3]==-1) {
      pars[,c(2,4)]<-NA
      errorflag<-1
    } else {
      errorflag<-0
    } # end if
  } else if(typ=='lm') {
    if(res$errorflag==-1) {
      pars[,c(2,4)]<-NA
      errorflag<-1
    } else {
      errorflag<-0
    } # end if
  } # end if
  ### Calculate signal values for each component
  CompSig<-apply(cbind(pars[,1],pars[,3]),MARGIN=1,
           function(x) if(typ=='cw') x[1]*exp(-x[2]*tim) else 
           x[1]*x[2]*(tim/max(tim))*exp(-x[2]*tim^2/2/max(tim)))
  CompSig<-round(cbind(res$predtval,CompSig),5L)
  colnames(CompSig)<-c('Fit.Signal', paste('Comp.',1:ncomp,sep='') )
  ### plot or not
  if(plot==TRUE) {
    ### set default xylim
    if(is.null(xylim)) {
      xylim<-if(typ=='cw') c(5,500) else c(1000,500)
    } # end if
    ### add a scatter plot (Time .VS. Signal)
    plot(tim,sig,main=samplename,
         xlab='Time (s)', ylab='OSL Counts',xlim=c(0,xylim[1]),ylim=c(0,xylim[2]),
         xaxs='i',yaxs='i',type='p',pch=21,cex=1,bg='grey',col='black',lty=1)
    ### set colors
    col<-c('black','blue')
    ### lines Time .VS. Fitted values
    lines(tim,CompSig[,1],lwd=lwd,col=col[1],lty=1)
    ### lines Time .VS. Component signal (1 to ncomp)
    for(i in 1:(ncomp-1)) {
        lines(tim,CompSig[,i+1],lwd=lwd,col=col[2],lty=i)
      } # end for
    ### lines Time .VS. Back.ground signal
    lines(tim,CompSig[,ncomp+1],lwd=lwd,col=col[1],lty=3)
    ### add a legend
    legend(ifelse(typ=='cw',"topright","topleft"),legend=c('Fitted.Curve',paste('Comp.',1:ncomp,sep='')),
           col=col[c(1,rep(2,ncomp-1),1)],pch=c(21,rep(NA,ncomp)),lty=c(1,1:(ncomp-1),3),
           yjust=2,ncol=1,cex=1,bty="o",lwd=lwd,pt.bg='grey')
  } # end if
  ### results for output
  out<-list('Comp.Signal'=CompSig,'pars'=pars,'value'=res$value,
            'errorflag'=errorflag)
  ### wirte fitted signal values to a file or not
  if(!is.null(outfile)) {
    write.csv(CompSig,file=paste(outfile,'.csv'))
  } # end if
  ### out is invisible
  invisible(out)
} # end function decomp
######################################## END FUNCTION decomp ####################################
########################################## FUNCTION sgcED #######################################
### ***********************************************************************************************
### Function sgcED is used to analyze equivalent doses using SGC method.
###
###     Author: Peng Jun, 2013.06.23, revised in 2013.07.26, revised in 2013.08.02
###
### References: Roberts,H.M. and Duller,G.A.T., 2004. Standardised growth curves for optical dating 
###             of sediment using multiple-grain aliquots. Radiation Measurements 38, pp. 241-252.
###
###             Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from 
###             single aliquot regenerative dose measurements. Ancient TL 25, pp. 15-24.
### ***********************************************************************************************
sgcED<-
function(Curvedata,Ltx,model=c('line','exp','line+exp'),
         nstart=100,upb=5,ErrorMethod=c('mc','sp'),MCiter=1000,
         plot=TRUE,samplename=NULL,outfile=NULL)  {
   UseMethod('sgcED')
}
### default method for function sgc.ED
sgcED.default<-
function(Curvedata,Ltx,model=c('line','exp','line+exp'),
         nstart=100,upb=5,ErrorMethod=c('mc','sp'),MCiter=1000,
         plot=TRUE,samplename=NULL,outfile=NULL)  {
  ### stop if not
  if(!is.data.frame(Curvedata))         stop('Error: Curvedata must be of type data.frame!')
  if(ncol(Curvedata)!=3L)               stop('Error: Curvedata must have three columns!')
  if(!is.numeric(
     unlist(unclass(Curvedata))))       stop('Error: elements in Curvedata must be all of type numeric!')
  if(!is.data.frame(Ltx))               stop('Error: Ltx must be a data.frame!')
  if(ncol(Ltx)!=2L)                     stop('Error: Ltx must contains two columns!')
  if(!is.numeric(
     unlist(unclass(Ltx))))             stop('Error: elements in Ltx must be all of type numeric!')
  if(any(Ltx<=0))                       stop('Error: Ltx must larger than 0!')
  if(any(Ltx[,1]>=max(Curvedata[,2])))  stop('Error: Ltx must not beyond maximum Lx/Tx in Curvedata!')
  if(!is.character(model))              stop('Error: model must be of type character!')
  if(length(model)>=4L)                 stop('Error: model must contain no more than three elements!')
  if(length(model)==1L) {
    if(!model %in% 
       c('line','exp','line+exp'))      stop('Error: model must be one of "line", "exp", "line+exp"!')
  } else if(length(model)>=2L) {
    if(!all(model %in% 
       c('line','exp','line+exp')))     stop('Error: incorrect model!')
  } # end if
  if(!is.numeric(nstart))               stop('Error: nstart must be of type numeric!')
  if(length(nstart)!=1L)                stop('Error: nstart must be an one-element vector!')
  if(abs(nstart-round(nstart))
     >=.Machine$double.eps^0.5)         stop('Error: nstart must be an integer!')
  if(nstart<=1L)                        stop('Error: nstart must larger than 1!')
  if(nstart>10000L)                     stop('Error: nstart is too large!')          
  if(!is.numeric(upb))                  stop('Error: upb must be of type numeric!')
  if(length(upb)!=1L)                   stop('Error: upb must be an one-element vector!')
  if(upb<=0)                            stop('Error: upb must larger than 0!')
  if(upb>1e4)                           stop('Error: upb is too large!')
  if(!is.character(ErrorMethod))        stop('Error: ErrorMethod must be of type character!')
  if(length(ErrorMethod)>=3L)           stop('Error: ErrorMethod must contain no more than two elements!')
  if(length(ErrorMethod)==1L) {
    if(!ErrorMethod %in%
       c('sp','mc'))                    stop('Error: ErrorMethod must be either "sp" or "mc"!')
  } else if(length(ErrorMethod==2L)) {
    if(!all(ErrorMethod %in% 
       c('sp','mc')))                   stop('Error: incorrect ErrorMethod!')
  } # end if
  if(!is.numeric(MCiter))               stop('Error: MCiter must be of type numeric!')
  if(length(MCiter)!=1L)                stop('Error: MCiter must be an one-element vector!')
  if(abs(MCiter-round(MCiter))
     >=.Machine$double.eps^0.5)         stop('Error: MCiter must be an integer!')
  if(MCiter<100L)                       stop('Error: MCiter is too small!')
  if(MCiter>10000L)                     stop('Error: MCiter is too large!')
  if(!is.logical(plot))                 stop('Error: plot must be either TRUE or FALSE!')
  if(!is.character(samplename) &&
     !is.null(samplename))              stop('Error: samplename must be NULL or of type character!')
  if(!is.null(outfile) &&
     !is.character(outfile))            stop('Error: outfile must be NULL or of type character!')
  ### set default model (linear)
  model<-model[1]
  ### set default ErrorMethod ('mc')
  ErrorMethod<-ErrorMethod[1]
  ### check if data is enough for model fitting
  if(model=='line' && nrow(Curvedata)<2L) {
    stop('Error: fitting a linear model needs at least two paired observations')
  } # end if
  if(model=='exp' && nrow(Curvedata)<3L) {
     stop('Error: fitting a exponential model needs at least three paired observations')
  } # end if
  if(model=='line+exp' && nrow(Curvedata)<4L) {
    stop('fitting a linear+exponential model needs at least four paired observations')
  } # end if
  ### parameters for subroutine 'sgcED.f90'
  Dose<-drop(Curvedata[,1])
  ltx<-cbind(drop(Curvedata[,2]),drop(Curvedata[,3]))
  ndose<-length(Dose)
  inltx<-cbind(drop(Ltx[,1]),drop(Ltx[,2]))
  ninltx<-nrow(Ltx)
  outDose<-matrix(0,nrow=ninltx,ncol=2L)
  npars<-if(model=='line') {
    2 } else if(model=='exp') {
    3 } else if(model=='line+exp') {
    4 } # end if
  pars<-parserrors<-vector(length=npars)
  predtval<-vector(length=ndose)
  value<- -99.0
  method<-ifelse(ErrorMethod=='sp',1,2)
  errorflag<-vector(length=2L)
  ### calculate equivalent doses
  res<-.Fortran('sgcED',as.double(Dose),as.double(ltx),as.integer(ndose),as.integer(ninltx),as.double(inltx),
        outDose=as.double(outDose),pars=as.double(pars),as.integer(npars),predtval=as.double(predtval),
        as.double(upb),as.integer(nstart),parserrors=as.double(parserrors),value=as.double(value),
        as.integer(method),as.integer(MCiter),errorflag=as.integer(errorflag))
  ### error checking
  if(res$errorflag[1]!=123) {
    stop('Error: fail in Dose-Response curve fitting, the fitting model might be inappropriate!')
  } # end if
  ### set LMpars
  LMpars<-cbind(res$pars,res$parserrors)
  colnames(LMpars)<-c('Estimate','Std.Error')
  rowname<-c('a','b','c','d')
  rownames(LMpars)<-rowname[1:npars]    
  if(res$errorflag[2]!=0) {
    LMpars[,2]<-NA
  } # end if
  ### set fit.value for output
  fit.value<-cbind(drop(Curvedata[,1]),drop(Curvedata[,2]),res$predtval)
  colnames(fit.value)<-c('ReDose','Lx/Tx','Fit.Lx/Tx')
  rownames(fit.value)<-paste('ReDose.',1:ndose,sep='')
  ### set ED values
  ED<-matrix(res$outDose,ncol=2L)
  rownames(ED)<-paste('NO.',1:ninltx,sep='')
  colnames(ED)<-c('ED','Std.ED')
  ### prepare results for output
  output<-list('LMpars'=LMpars,'residual'=res$value,
               'fit.value'=fit.value,'ED'=ED)
  ### plot or not
  if(plot==TRUE) {
    ### set limitations for x and y axies
    Xlim<-max(Curvedata[,1],ED[,1])
    Ylim<-max(Curvedata[,2],Ltx[,1])
    ### plot a dose-response curve
    plot(NA,NA,main=samplename,
         xlab='ReDose (Gy)',ylab='Lx/Tx',
         xlim=c(0,Xlim*1.05),ylim=c(0,Ylim*1.05),
         xaxs='i',yaxs='i',lab=c(7,7,9) )
    ### add ReDose as points
    points(Curvedata[,1:2],pch=1,cex=3)
    ### points calculate Equivalent Dose .VS. Ltx
    points(ED[,1],Ltx[,1],pch=23,cex=3,bg='grey')
    ### add error bars to Lx/Tx
    if(all(Curvedata[,3]>=1e-3)) {
      arr<-suppressWarnings(try(arrows(drop(Curvedata[,1]),drop(Curvedata[,2])-drop(Curvedata[,3])/2,drop(Curvedata[,1]),
           drop(Curvedata[,2])+drop(Curvedata[,3])/2,code=3,lwd=2.5,angle=90,length=0.05,col="black"),silent=TRUE))
    } # end if
    ### add error bars to calculated ED
    if(all(ED[,2]>=1e-3)) {
      arr<-suppressWarnings(try(arrows(ED[,1]-ED[,2]/2,Ltx[,1],ED[,1]+ED[,2]/2,
           Ltx[,1],code=3,lwd=2.5,angle=90,length=0.05,col="black"),silent=TRUE))
    } # end if
    ### reset model
    model<-model[1]
    ### add a fitting curve
    x<-NULL
    if(model=='line') {
      curve(LMpars[1,1]*x+LMpars[2,1],type='l',
            add=TRUE,lw=2,from=0,to=Xlim*1.05)
    } else if(model=='exp') {
      curve(LMpars[1,1]*(1-exp(-LMpars[2,1]*x))+LMpars[3,1],type='l',
            add=TRUE,lw=2,from=0,to=Xlim*1.05)
    } else if(model=='line+exp')  {
      curve(LMpars[1,1]*(1-exp(-LMpars[2,1]*x))+LMpars[3,1]*x+LMpars[4,1],type='l',
            add=TRUE,lw=2,from=0,to=Xlim*1.05)
    } # end if
    ### add dash lines 
    for(i in 1: ninltx) {
      lines(c(0,ED[i,1],ED[i,1]),c(Ltx[i,1],Ltx[i,1],0),lty='dashed',lwd=0.5)
    } # end for
  } # end if
  ### write equivalent doses out or not
  if(!is.null(outfile)) {
    write.csv(ED,file=paste(outfile,'.csv'))
  } # end if
  return(output)
} # end function sgcED   
######################################## END FUNCTION sgcED ###############################################
