# Wrapped function RadialPlotter calls 
# function RadialPlotter.default
RadialPlotter<-
function(EDdata,ncomp=0,addsigma=0,
         maxiter=500,maxcomp=9,
         eps=.Machine$double.eps^0.5,plot=TRUE,
         zscale=NULL,samplename=NULL)  {
  UseMethod('RadialPlotter')
}
#
############################################### RadialPlotter.default ####################################################
#
# set default method for function RadialPlotter
RadialPlotter.default <-
function(EDdata,ncomp=0,addsigma=0,
         maxiter=500,maxcomp=9,
         eps=.Machine$double.eps^0.5,plot=TRUE,
         zscale=NULL,samplename=NULL)  {
  if(!is.data.frame(EDdata))  stop('EDdata must be type of data.frame!')
  if(any(EDdata[,1]<=0))      stop('Equivalent Doses must be larger than 0!')
  if(!ncomp%in%(-2:maxcomp))  stop('ncomp must be an integer ranging from -2 to maxcomp!')
  if(ncomp==0&&maxcomp<=1)    stop("maxcomp must non't be smaller than 2 if ncomp is 0!") 
  if(addsigma<0)  stop('addsigma must not be smaller than 0!')
  if(maxcomp<0)   stop('maxcomp must not be smaller than 0!')
  n<-nrow(EDdata)
  if(maxcomp>n)   stop('maxcomp must not exceed the length of ED!')
  ed<-EDdata[,1]
  error<-EDdata[,2]
  #
  # function for radial plot drawing,
  # the code is revised from Rex Galbraith 
  RadialPlot<-function(Data,Pars,zscale,samplename)  {
    z.i<-log(Data[,1])
    se.i<-Data[,2]/Data[,1]
    # if Pars=NULL, lines will not be plotted out
    if(!is.null(Pars))  {
      Pars<-log(Pars)
    }
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
      }
    } else  {
      mkest<-zscale
    } 
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
      12   }  
    par(mfrow=c(1,1),oma=c(2,0,0,0),mar=c(5,4,4,4.5),xpd=TRUE,las=1,cex=1)
    plot(NA,NA,xlim=c(0,max(x.i)),ylim=c(-yaxis.max,yaxis.max),xlab='',ylab='',
       xaxt='n',yaxt='n',main=samplename,bty='n',typ='n',xaxs='i',yaxs='i')      
    lenpar<-length(Pars) 
    maxx<-which.max(circ.x)
    # if Pars!=NULL, now starting drawing lines out
    if(!is.null(Pars))  {
      for(i in 1:lenpar)  {
        if(Pars[i]-z.0>=0)  {
          Ci<-approx(circ.y[maxx:201],circ.x[maxx:201],ties="ordered",
          xout=(Pars[i]-z.0)*max(x.i),rule=2,method="constant",f=0.5)$y
        }  else  {
          Ci<-approx(circ.y[1:(maxx-1)],circ.x[1:(maxx-1)],ties="ordered",
          xout=(Pars[i]-z.0)*max(x.i),rule=2,method="constant",f=0.5)$y
        }
       lines(c(0,Ci),c(0,(Pars[i]-z.0)*Ci),lty=1,col='black',lwd=1.5)
      }   
    }
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
  } # end function RadialPlot
  #
  tol<-.Machine$double.eps^0.3    
  errorflag<-0
  maxlik<-BIC<-0
  #
  # for CAM and FMM analysis
  if(ncomp%in%(0:maxcomp))  {
    if (ncomp==0)  {
      goodcomp<-0
      BILI<-matrix(0,nrow=maxcomp-1,ncol=2L)
      # call Fortran subroutine 'FineComp' to pick out the appropriate 
      # number of component automatically
      Results<-.Fortran('FineComp',as.double(ed),as.double(error),as.integer(n),
                goodcomp=as.integer(goodcomp),as.double(addsigma),as.integer(maxiter),       
                as.double(eps),as.integer(maxcomp),BILI=as.double(BILI),pacakge='RadialPlotter')
      BILI<-Results$BILI
      dim(BILI)<-c(maxcomp-1,2L)
      colnames(BILI)<-c('BIC','Maxlik')
      rownames(BILI)<- paste(rep("k=", maxcomp-1),c(2:maxcomp),sep = "") 
      ncomp<-Results$goodcomp
      loopcomp<-1
    }  
    spars<-pars<-matrix(0,nrow=2L,ncol=ncomp)
    # call Fortran subroutine 'FineED1' to calculate parameters in 
    # central age model or finite mixture age model
    Results<-.Fortran('FineED1',as.double(ed),as.double(error),as.integer(n),
                      as.integer(ncomp),as.double(addsigma),
                      pars=as.double(pars),spars=as.double(spars),maxlik=as.double(maxlik),
                      BIC=as.double(BIC),as.integer(maxiter),as.double(eps),
                      as.double(tol),errorflag=as.integer(errorflag),pacakge='RadialPlotter')  
    ParsAndErrors<-cbind(matrix(Results$pars,byrow=TRUE,ncol=2),matrix(Results$spars,byrow=TRUE,ncol=2))
    if(ncomp==1)  { # for CAM analysis
      ParsAndErrors<-matrix(ParsAndErrors[,c(1,3,2,4)])
      rownames(ParsAndErrors)<-c('Sigma','sSigma','ED','sED')
      colnames(ParsAndErrors)<-'CAM'
      if(plot)  {
        RadialPlot(EDdata,ParsAndErrors[3],zscale,samplename)
      }
    }  else {      # for FMM analysis
      ParsAndErrors<-ParsAndErrors[,c(1,3,2,4)][order(ParsAndErrors[,2]),]
      if(plot)  {
        RadialPlot(EDdata,ParsAndErrors[,3],zscale,samplename)
      }
      colnames(ParsAndErrors)<-c('P','sP','ED','sED') 
      rownames(ParsAndErrors) <- paste(rep("comp", ncomp), c(1:ncomp), sep = "") 
    } 
  }  else  { 
    # for MAM analysis 
    pars<-spars<-vector(length=2-ncomp)
    # call Fortran subroutine 'MAM' to fit the minimum age models
    Results<-.Fortran('MAM',as.double(ed),as.double(error),as.integer(n),pars=as.double(pars),
                     spars=as.double(spars), maxlik=as.double(maxlik), BIC=as.double(BIC), 
                     as.integer(2-ncomp),as.double(addsigma), as.integer(maxiter),
                     as.double(tol),errorflag=as.integer(errorflag),pacakge='RadialPlotter')
    ParsAndErrors<-cbind(Results$pars,Results$spars)
    colnames(ParsAndErrors)<-c('pars','error')
    if(ncomp==-1)  {
      rownames(ParsAndErrors)<-c('P','gama','sigma')
    }  else  {
      rownames(ParsAndErrors)<-c('P','gama','mu','sigma')
    }
    if(plot)  {
      if(any(Results$spars==0)==FALSE)  {
        RadialPlot(EDdata,Results$pars[2],zscale=zscale,samplename=samplename)
      }  else {
        RadialPlot(EDdata,Pars=NULL,zscale=zscale,samplename=samplename)
      }
    }
  }
  # at this point the age model analysis have been terminated
  # now calculate the commom age model based equivalent dose
  z.i<-log(ed)
  se.i<-sqrt( (error/ed)^2 + addsigma^2 )
  commonED<-sum(z.i/se.i^2)/sum(1/se.i^2)
  scommonED <- 1/sqrt(sum(1/se.i^2))
  # the returned list
  out<-list(errorflag=Results$errorflag,commonED=exp(commonED)*c(1,scommonED),
            pars=ParsAndErrors,BIC=Results$BIC,maxlik=Results$maxlik)
  invisible(out)
  # now output the result to the terminal screen
  cat('\n')
  cat('================Results for RadialPlotter================','\n\n')
  cat('Error message:',Results$errorflag,'\n\n')
  cat('Common equivalent dose:',round(exp(commonED),4L),'+-',round(exp(commonED)*scommonED,4L),'\n\n')
  if(ncomp==1)  {
    cat('Overdispersion and Equivalent Dose are:','\n\n')
  }  else if(ncomp%in%c(0,2:maxcomp))  {
    if(exists('loopcomp'))  {
      cat('The best component number is:',ncomp,'\n\n')
    }
    cat('Parameters and standard errors are:','\n')
  }  else  {
    if(any(Results$spars==0)==FALSE)  {
      if(Results$errorflag==1) {
        cat('Warning: at least one parameter is near to the boundary!','\n\n')
      }
      cat('Parameters and standard errors are:','\n')
    }  else  {
      cat('Error: parameters and standard errors can not be estimated!','\n\n')
    }
  }
  if(!ncomp%in%(-2:-1) || any(Results$spars==0)==FALSE)  {
    cat('----------------------------------','\n')
    print(round(as.data.frame(ParsAndErrors),4L))
    cat('----------------------------------','\n\n')
    cat('                      BIC value:',round(Results$BIC,4L),'\n\n')
    cat('Maximum logged likelihood value:',round(Results$maxlik,4L),'\n\n')
    if(exists('loopcomp'))  {
      cat('BICs and Maximum logged likelihood values are:','\n')
      cat('------------------------','\n')
      print(round(as.data.frame(BILI),4L))
      cat('------------------------','\n\n')
    }
  }
  cat('=========================The End=========================','\n\n')
  }  # end function RadialPlotter.default
  #
  #
######################################################## The End ###################################################################
