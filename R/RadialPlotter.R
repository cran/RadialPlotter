###
# Wrapped function RadialPlotter calls 
# function RadialPlotter.default
################################################
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
  if(any(EDdata[,1]<=0))      stop('Equivalent must be larger than 0!')
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
################################################### End FUNCTION RadialPlotter #########################################################
###
######
### Author: Peng Jun, 2013.05.28
###
### References: Roberts,H.M. and Duller,G.A.T., 2004. Standardised growth curves for optical dating of sediment
###             using multiple-grain aliquots. Radiation Measurements 38, pp. 241-252.
###
###             Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from single aliquot
###             regenerative dose measurements. Ancient TL 25, pp. 15-24.
###
################################################### FUNCTION SGCED ######################################################## 
SGCED<-function(CurveData,
                Ltx,
                model='exp',
                iniPars=NULL,
                ErrorMethod='mc',
                MCiter=1e4,
                plot=TRUE,
                dashline=TRUE)   {
   UseMethod('SGCED')
} # end function SGCED
###
###########################################
###
### default method for function SGCED
SGCED.default<-function(CurveData,
                        Ltx,
                        model='exp',
                        iniPars=NULL,
                        ErrorMethod='mc',
                        MCiter=1e4,
                        plot=TRUE,
                        dashline=TRUE)        {
  ###
  ### argument checking for function SGCED
  stopifnot(nrow(CurveData)>=2,
            !is.null(ncol(CurveData)),
            ncol(CurveData)==2,
            all(CurveData>=0),
            all(Ltx>0),
            MCiter>=1e2,
            MCiter<=1e5,
            is.numeric(iniPars) || is.null(iniPars),
            model%in%c('line','exp','line+exp'),
            ErrorMethod%in%c('sp','mc'),
            is.logical(plot),
            is.logical(dashline)
            )
  ###
  ### 
  if(is.null(iniPars)) {
    # initilize parameters if initial 
    # parameters have not been provided
    if(model=='line') {
      iniPars<-c(1,1)
    } else if(model=='exp') {
      iniPars<-c(5.0,0.05,0.05)
    } else {
      iniPars<-c(5.0,0.05,0.05,0.05)
    } # end if
  } else {
    # check if iniPars have been provided appropriately
    if(model=='line' && length(iniPars)!=2) {
      stop('Linear Model Needs Two Initial Parameters')
    } # end if
    if(model=='exp' && length(iniPars)!=3) {
      stop('Exponential Model Needs Three Initial Parameters')
    } # end if
    if(model=='line+exp' && length(iniPars)!=4) {
      stop('Linear+Exponential Model Needs Four Initial Parameters')
    } # end if
  } # end if
  ###
  ###
  # check if data is enough for model fitting
  ndose<-nrow(CurveData)                # ndose
  if(model=='line' && ndose<2) {
    stop('Fitting a Linear Model Meeds at Least Two Paired Observations')
  } # end if
  if(model=='exp' && ndose<3) {
     stop('Fitting a Exponential Model Needs at Least Three Paired Observations')
  } # end if
  if(model=='line+exp' && ndose<4) {
    stop('Fitting a Linear+Exponential Model Needs at Least Four Paired Observations')
  } # end if
  ###
  ### ****************************************************************************
  ### specify parameters that will be used 
  ### in Fortran subroutine 'sgc.f90'
  dose<-CurveData[,1]                    # dose
  ltx<-CurveData[,2]                     # ltx
  ###
  ### reshape inltx
  if(is.null(ncol(Ltx))) {
    inltx<-cbind(Ltx,0)
  } else if(ncol(Ltx)==2) {              # inltx
    inltx<-cbind(Ltx[,1],Ltx[,2])
  } else {
    stop('Incorrect Format of Ltx')
  } # end if
  nout<-nrow(inltx)                      # nout
  outDose<-matrix(0.0,nrow=nout,ncol=2)  # outDose                  
  pars<-iniPars                          # pars
  predtval<-vector(length=ndose)         # predtval
  npars<-length(iniPars)                 # npars
  parserrors<-vector(length=npars)       # parserrors 
  value<--99.0                           # value
  errorflag<-vector(length=2)            # errorflag
  ###
  ### set method
  if(ErrorMethod=='sp')  {
    method<-1
  } else {                               # method
    method<-2
  } # end if
  ### ****************************************************************************
  ###
  ### calling Fortran subroutine 'sgc.f90'
  res<-.Fortran('sgc',as.double(dose),as.double(ltx),as.integer(ndose),as.double(inltx),
                      outDose=as.double(outDose),as.integer(nout),pars=as.double(pars),
                      predtval=as.double(predtval),parserrors=as.double(parserrors),
                      as.integer(npars),value=as.double(value),errorflag=as.integer(errorflag),
                      as.integer(method),as.integer(MCiter),package='RadialPlotter'
               )
  ###
  ### redefine LMpars for output
  LMpars<-cbind(res$pars,res$parserrors)
  if(res$errorflag[2]!=0) {
    LMpars[,2]<-NA
  } # end if
  colnames(LMpars)<-c('Estimate','Std.Error')
  rowname<-c('a','b','c','d')
  rownames(LMpars)<-rowname[1:npars]
  ###
  ### redefine fit.value for output
  fit.value<-cbind(CurveData,res$predtval)
  colnames(fit.value)<-c('ReDose','Lx/Tx','Fit.Lx/Tx')
  rownames(fit.value)<-paste('ReDose.',1:ndose,sep='')
  ###
  ### redefine Dose for output
  Dose<-matrix(res$outDose,ncol=2)
  if( is.null(ncol(Ltx))  ) {
    Dose[,2]<-NA
  } # end if
  colnames(Dose)<-c('ED','Std.ED')
  rownames(Dose)<-paste('ED',1:nout,sep='')
  ###
  ###
  ### define error message for output
  error<-vector(length=4)
  ### 1) check error message in Levenberg-Marquadt fitting
  if(res$errorflag[1]==123) {
    error[1]<-'LM_Ok'
  } else {
    error[1]<-'LM_Fail'
  } # end if
  ### 2) check if characteristic parameters' standard errors can be estimated
  if(res$errorflag[2]==0)  {
    error[2]<-'LM.Std_Ok'
  } else {
    error[2]<-'LM.Std_Fail'
  } # end if 
  ### 3) check if Dose can be calculated
  if( any(!is.finite(Dose[,1]) ) || 
      any( abs(na.omit(Dose[,1])+99.0)<=.Machine$double.eps^0.1 )  ) {
    error[3]<-'ED_Fail'
  } else {
    error[3]<-'ED_Ok'
  } # end if
  ## 4) check if Doses' standard errors can be estimated
  if( any(!is.finite(Dose[,2]))  )  {
    error[4]<-'ED.Std_Fail'
  } else {
    error[4]<-'ED.Std_Ok'
  } # end if
  ###
  ### remove NA, Inf from inltx and Dose
  Findex<-which(is.finite(Dose[,1]) )              
  inltx<-inltx[Findex,,drop=FALSE]
  Dose<-Dose[Findex,,drop=FALSE]
  ### remove -99.0 from inltx and Dose
  Findex<-which( abs(Dose[,1]+99.0)>.Machine$double.eps^0.1 )
  inltx<-inltx[Findex,,drop=FALSE]
  Dose<-Dose[Findex,,drop=FALSE]
  ### set those minus standard errors of Dose to be NA
  NAS<-which(Dose[,2]<0)
  Dose[NAS,2]<-NA
  ### check if any number exist in Dose 
  if(nrow(Dose)<1) {
    Dose<-NULL
  } # end if
  ###
  ### diagnose if the fitting is successful
  if(error[1]=='LM_Fail' ||
     error[2]=='LM.Std_Fail' ) {
    cat('Warning: Bad Estimation: Inappropriate Fitting Model or Improper Initial Guess Values!','\n\n')
  } # end if
  ###
  ###
  ### check how many ED values can not be calculated
  if(error[3]=='ED_Fail') {
    cat(paste('Warning: Leaving',ifelse(is.null(Dose),'all',nout-nrow(Dose)),
              'Equivalent Dose Values Uncalculated!',sep=' '),'\n\n')
  } # end if
  ### prepare results for output
  output<-list('errormessage'=as.matrix(error),
               'LMpars'=LMpars,
               'residual'=res$value,
               'fit.value'=fit.value,
               'Dose'=Dose )
  ###
  ### plot or not
  if(plot)  {
    ### perform the plot only if Dose data is available
    if(!is.null(Dose)) {
      ###
      ### plot paired observations
      Xlim<-max(CurveData[,1],Dose[,1])
      Ylim<-max(CurveData[,2],inltx[,1])
      plot(CurveData,
           main='Dose-Response Curve',
           xlab='ReDose (Gy)',
           ylab='Lx/Tx',
           xlim=c(0,Xlim*1.05),
           ylim=c(0,Ylim*1.05),
           pch=18,
           cex=3,
           xaxs='i',
           yaxs='i',
           lab=c(7,7,9) )
      ###
      ### points out calculate Dose values corresponded to Ltx
      points(Dose[,1],inltx[,1],pch=5,cex=2,col='blue')
      ###
      ### add a fitting curve to the plot
      x<-NULL
      if(npars==2) {
        curve(LMpars[1,1]*x+LMpars[2,1],type='l',add=TRUE,lw=2,from=0,to=Xlim*1.05)
      } else if(npars==3) {
        curve(LMpars[1,1]*(1-exp(-LMpars[2,1]*x))+LMpars[3,1],type='l',add=TRUE,lw=2,from=0,to=Xlim*1.05)
      } else {
        curve(LMpars[1,1]*(1-exp(-LMpars[2,1]*x))+LMpars[3,1]*x+LMpars[4,1],type='l',add=TRUE,lw=2,from=0,to=Xlim*1.05)
      } # end if
      ###
      ### add a legend
      legend("bottomright", 
             legend=paste('Fitting Model:',model,sep=' '),
             col='black',
             yjust=2,
             ncol=1,
             bty="n")
      ###
      ### and dash line for each calculated paired value 
      if(dashline) {
        for(i in 1:nrow(Dose) ) {
          lines(c(0,Dose[i,1],Dose[i,1]),c(inltx[i,1],inltx[i,1],0),lty='dashed',lwd=1)
        } # end for
      } # end if
      ###
      ### add error bars only if standard errors for Doses are available
      Findex<-which(!is.na(Dose[,2]))    # here Findex has been reset
      arr<-suppressWarnings(arrows(Dose[Findex, 1]-Dose[Findex, 2]/2,inltx[Findex,1],Dose[Findex, 1]+Dose[Findex, 2]/2,
                            inltx[Findex,1],code=3,lwd=2.5,angle=90,length=0.05,col="black"))
    } else {
      cat('Warning: Failing in Plot Drawing, No Available Equivalent Dose!','\n\n')
    } # end if
  } # end if
  ### 
  ###output the results    
  return(output)
  ###
} # end function SGCED.default
###
################################################ END FUNCTION SGCED ######################################################## 
###
###
###
### Author: Peng Jun, 2013.06.05
###
####################################
###
### References : Bluszcz, A., 1996. Exponential function fitting to TL growth data and similar
###              applications. Geochronometria 13, 135–141.
###
###              Bluszcz, A., Adamiec, G., 2006. Application of differential evolution to fitting
###              OSL decay curves. Radiation Measurements 41, 886-891.
###   
###              Jain, M., Murray, A.S., Boetter-Jensen, L., 2003. Characterisation of blue-light stimulated 
###              luminescence components in different quartz samples: implications for dose measurement. Radiation
###              Measurements, 37 (4-5), pp. 441-449.
### 
######################################### FUNCTION DECOMP ###############################################
###
decomp<-function(SigData,
                 ncomp=3,
                 LEDpower=60,
                 LEDwavelength=470,
                 plot=TRUE,
                 xyLim=c(5,500),
                 typ='CW',
                 control.args=list() )  {
  ###
  UseMethod('decomp')
  ###
} # end function decomp
###
#####################################
###
### default method for decomp
###
decomp.default<-function(SigData,
                         ncomp=3,
                         LEDpower=60,
                         LEDwavelength=470,
                         plot=TRUE,
                         xyLim=c(5,500),
                         typ='CW',
                         control.args=list() )  {
  ###
  stopifnot(is.numeric(ncol(SigData)),
            ncol(SigData)==2,
            all(SigData>0.0),
            ncomp%in%(1:7),
            is.numeric(LEDpower),
            is.numeric(LEDwavelength),
            is.numeric(xyLim),
            all(xyLim>0.0),
            is.logical(plot),
            typ=='CW',
            class(control.args)=='list',
            all(names(control.args) %in% list('factor','f','cr','maxiter','tol') )
            )
  ###
  tim<-SigData[,1]
  sig<-SigData[,2]  
  ntim<-length(tim)
  pars<-Stdpars<-vector(length=2*ncomp)
  value<--99.0
  predtval<-vector(length=ntim)
  errorflag<-vector(length=3)
  ###
  ###
  ### default arguments for differential evolution
  ### ****************************************************
  args<-list(factor=10,
             f=0.5,
             cr=0.99,
             maxiter=1e3,
             tol=0.1)   
  args[names(control.args)]<-control.args
  factor<-args$factor
  f<-args$f
  cr<-args$cr
  maxiter<-args$maxiter
  tol<-args$tol
  if(factor<3) {
    stop('control.args: factor is too small!')
  } # end if
  if(f<=0.0 || f>1.2) {
    stop('control.args: incorrect f!')
  } # end if
  if(cr<=0.0 || cr>1.0) {
    stop('control.args: incorrect cr!')
  } # end if
  if(maxiter<10) {
    stop('control.args: maxiter is too small!')
  } # end if
  if(maxiter>1e4) {
    stop('control.args: maxiter is too large!')
  } # end if
  if(tol<0.0) {
    stop('control.args: incorrect tol!')
  } # end if
  ### ****************************************************
  ###
  ###
  res<-.Fortran('decomp',as.integer(ncomp),as.double(tim),as.double(sig),
                         as.integer(ntim),pars=as.double(pars),Stdpars=as.double(Stdpars),
                         value=as.double(value),predtval=as.double(predtval),
                         as.integer(factor),as.double(f),as.double(cr),as.integer(maxiter),
                         as.double(tol),errorflag=as.integer(errorflag),package='RadialPlotter' 
                )
  ###
  ### Error checking
  if(res$errorflag[1]==1 && 
     res$errorflag[3]==-1)  {
     stop(paste('Signal can not be decomposed to',ncomp,'Components!'))
  } # end if
  ###
  ###
  ### decide Photoionisation cross-section (cm2)
  h<-6.62606957e-34
  ny<-299792458/(LEDwavelength/10^9)
  E<-h*ny
  LEDpower<-LEDpower/1000
  ###
  ###
  ### reshpae pars for output
  pars<-cbind(res$pars[1:ncomp],
              res$Stdpars[1:ncomp], 
              res$pars[-(1:ncomp)],
              res$Stdpars[-(1:ncomp)],
              res$pars[-(1:ncomp)]/LEDpower*E)
  pars<-pars[order(pars[,3],decreasing=TRUE),,drop=FALSE]
  colnames(pars)<-c('Ithn','Std.Ithn','Lamda','Std.Lamda','Pcs')
  rownames(pars)<-paste('Comp.',1:ncomp,sep='')
  ###
  ###
  ### reset errorflag
  if(res$errorflag[3]==-1) {
    pars[,2]<-NA
    pars[,4]<-NA 
    errorflag<-1
  } else {
    errorflag<-0
  } # end if
  ###
  ###
  ### Calculate signal values for each component
  CompSig<-apply(cbind(pars[,1],pars[,3]),
                 MARGIN=1,
                 function(x) x[1]*exp(-x[2]*tim) )
  CompSig<-round(cbind(res$predtval,CompSig),3L)
  colnames(CompSig)<-c('Fit.Signal', paste('Comp.',1:ncomp,sep='') )
  ###
  ###
  ### plot or not
  if(plot) {
    plot(tim,
         sig,
         main=paste('Ncomp=',ncomp,sep=''),
         xlab='Time (s)', 
         ylab='OSL Counts',
         xlim=c(0,xyLim[1]),
         ylim=c(0,xyLim[2]),
         xaxs='i',
         yaxs='i',
         type='o',
         pch=20,
         cex=2)
    pch<-c(0,2,5,7,9,10)
    lines(tim,CompSig[,1],col='blue',pch=16,type='o',cex=1,lwd=2)
    lines(tim,CompSig[,ncomp+1], pch=1,cex=1,type='o',lwd=1)
    for(i in 2:ncomp ) {
      lines(tim,CompSig[,i],pch=pch[i-1],lwd=1,cex=1,type='o')
    } # end for
    legend("topright",
           legend=c('Obeservations','Fitted.Value',paste('Comp.',1:ncomp,sep='')),
           pch=c(20,16,pch[1:(ncomp-1)],1),
           col=c('black','blue',rep('black',ncomp)),
           yjust=2,
           ncol=1,
           cex=1,
           bty="n",
           lwd=1)
  } # end if
  ###     
  ###
  ### output
  list('Comp.Signal'=CompSig, 
       'pars'=pars, 
       'value'=res$value,
       'errorflag'=errorflag)
  ###
} # end function decomp
###
###
######################################### END FUNCTION DECOMP ###############################################
###
