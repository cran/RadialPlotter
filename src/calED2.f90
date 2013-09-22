subroutine calED2(Dose,ltx,ndose,inltx,outDose,pars,npars,predtval,upb,&
                  nstart,parserrors,value,mcED,method,motoiter,errorflag)
!-------------------------------------------------------------------------------------------------------------------------------
! Fitting a dose-response curve and calculate equivalent dose by interpolation (origin version)
!
! ndose,               input:: integer, length of reDoses (Redose1, Redose2,...) or Ltx(Lx1/Tx1, Lx2/Tx2,...)
! nstart,              input:: integer, number of random trials
! upb,                 input:: real value, upper boundary for b value, b is generated in (0, upb)
! npars,               input:: integer, model used for fitting the dose-response curve:
!                              if npars=1, a linear model of the form: y=ax will be fitted, where ndat>=1
!                              if npars=2, a Exponential model of the form: y=a(1-exp(-bx))  will be fitted, where ndat>=2
!                              if npars=3, a linear+Exponential model of the form: y=a(1-exp(-bx))+cx will be fitted, where ndat>=3
! dose(ndose),         input:: real values, Redose data used to build a dose-response curve
! ltx(ndose,2),        input:: real values, standardlized OSL signal data used to build a dose-response curve
! inltx(2),            input:: real values, standardlized OSL signal values that used for calculating EDs and EDs' standard errors
! outDose(2),         output:: real values, calculated EDs and standard errors that correlated to standardlized OSL signal
! pars(npars),        output:: real values, calculated characteristical parameters of a dose-response curve
! parserrors(npars),  output:: real values, parameters' standard errors
! predtval(ndose),    output:: real values, fitted values that correspond to ltx
! value,              output:: real value, sum of the square of residual
! mcED(motoiter),     output:: real values, EDs generated using monte carlo simulation
! motoiter,            input:: integer, numbers of Moto Carlo simulations that using to estimate EDs' standard errors
! method,              input:: integer, method for calculating dose's standard error, 1 for simple method; 2 for Moto Carlo method
! errorflag(2),       output:: integer, error message generated when calling subroutine SGC:
!                              1) error appears when calling subroutine lmder1, errorflag(1) will be one in (0,4,5,6,7,8,9,10), 
!                                 else errorflag(1)=123; error appears when calling linearfit, errorflag(1) will be 0, else 123
!                              2) error appears when attempt to calculate parameters' standard errors, errorflag(2)=1, else
!                                 errorflag(2) will be 0
!
! Author:: Peng Jun, 2013.09.20, revised in 2013.09.21
!
! Dependence:: subroutine inipars2; subroutine lmfit2; subroutine interpolate2; subroutine linearfit; subroutine r8vec_normal
!
!              Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from single aliquot
!              regenerative dose measurements. Ancient TL 25, pp. 15-24.
!----------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::ndose
  integer(kind=4),intent(in)::npars
  integer(kind=4),intent(in)::motoiter
  integer(kind=4),intent(in)::method
  integer(kind=4),intent(in)::nstart
  real   (kind=8),intent(in)::upb
  real   (kind=8),intent(out)::value
  integer(kind=4),intent(out)::errorflag(2)
  real   (kind=8),dimension(ndose),intent(in)::Dose
  real   (kind=8),dimension(ndose,2),intent(in)::ltx
  real   (kind=8),dimension(2),intent(in)::inltx
  real   (kind=8),dimension(2),intent(out)::outDose
  real   (kind=8),dimension(npars),intent(out)::pars
  real   (kind=8),dimension(npars),intent(out)::parserrors
  real   (kind=8),dimension(ndose),intent(out)::predtval
  real   (kind=8),dimension(motoiter),intent(out)::mcED
  ! local variables
  real   (kind=8),parameter::lmtol=1.0D-07                 ! toleracne for stopping iterations for subroutine lmfit2
  real   (kind=8)::minvalue                                ! minimized value of subroutine interpolate2
  real   (kind=8)::maxDose                                 ! maximum Redose value
  real   (kind=8)::averageErr                              ! average value of the sum of the squared residual errors
  real   (kind=8),dimension(2)::lowup                      ! low and up boundary dose values
  real   (kind=8)::lowltx,upltx                            ! low and up boundary ltx value
  integer(kind=4)::seed                                    ! seed for normal distribution random number generation
  real   (kind=8),dimension(ndose)::ranltx                 ! random Lx/Tx
  real   (kind=8),dimension(npars)::Ranpars,Ranparserrors  ! characteristic parameters of simulated dose-response curve 
  real   (kind=8),dimension(ndose)::Ranpredtval            ! fitted values correspond to ranltx
  real   (kind=8)::ranvalue                                ! sum of the square of residual generaged with random simulation
  integer(kind=4),dimension(2)::Ranerrorflag               ! errorflag generated during random fitting
  integer(kind=4)::i,j,mccount                             ! iterative count
  real   (kind=8)::sumdose,sumdose2,mcDose                 ! cumulative sum of ED, sum of the square of ED, mc ED of
  !                                                        ! Moto Carlo ED values
  ! variables for subroutine inipars2
  real   (kind=8),dimension(1):: rand
  integer(kind=4),dimension(2)::info
  real   (kind=8),dimension(npars-1)::outpars
  real   (kind=8),dimension(3)::cpars
  real   (kind=8),dimension(npars)::lmpars,lmparserrors
  real   (kind=8),dimension(ndose)::lmpredtval 
  integer(kind=4),dimension(2)::lmerrorflag
  real   (kind=8)::bvalue,lmvalue,itervalue
  ! default returned values if error appears
  outDose=-99.0D+00
  mcED=-99.0D+00
  value=-99.0D+00
  pars=-99.0D+00
  parserrors=-99.0D+00
  predtval=-99.0D+00
  ! default returned errorflag if error appears
  errorflag(1)=0
  errorflag(2)=1
  ! initialize parameters
  if(npars==1) then  ! for linear model
    ! call linearfit to obtain the characterized parameters of a dose-response
    ! curve, return lmpars, lmparserrors, lmpredtval, lmvalue and lmerrorflag
    call linearfit(Dose,ltx(:,1),ndose,lmpars,npars,&
                   .true.,lmparserrors,lmpredtval,&
                   lmvalue,lmerrorflag) 
    pars=lmpars
    parserrors=lmparserrors
    predtval=lmpredtval
    value=lmvalue
    errorflag=lmerrorflag
  else if(npars==2 .or. npars==3) then ! for expentional or linear plus expentional model
    cpars=0.0D+00
    itervalue=1.0D+30
    call random_seed()
    loop: do i=1, nstart
      call random_number(rand)
      bvalue=upb*rand(1)
      call inipars2(bvalue,npars+1,ndose,dose,&
                    ltx(:,1),outpars,lmtol,info)
      ! error checking of subroutine inipars2
      if(info(1)==0 .and. info(2)==0) then
        if(npars==2)  then
          cpars(1)=outpars(1)
          cpars(2)=bvalue
        else if(npars==3) then
          cpars(1)=outpars(1)
          cpars(2)=bvalue
          cpars(3)=outpars(2)
        end if
        ! pass cpars to lmpars
        lmpars=cpars(1:npars)
      end if
      if(info(1)/=0 .or. info(2)/=0) cycle loop
      ! call lmfit2 to get the characterized parameters for Dose-Response 
      ! curve, return pars, parserrors,predtval,value, errorflag
      call lmfit2(Dose,ltx(:,1),ndose,lmpars,npars,&
                  .true.,lmparserrors,lmpredtval,&
                  lmvalue,lmtol,lmerrorflag) 
      ! check if any error appears when calling lmfit2
      ! store improved results only
      if(lmerrorflag(1)==123  .and. lmvalue<itervalue) then
        itervalue=lmvalue
        pars=lmpars
        parserrors=lmparserrors
        predtval=lmpredtval
        value=lmvalue
        errorflag=lmerrorflag
      end if
    end do loop
  end if
  ! error checking, if error appears in linearfit or lmfit2, return
  if(errorflag(1)/=123) return
  ! calculate equivalent dose using inltx(1), store in outDose(1)
  maxDose=maxval(Dose)
  call interpolate2(outDose(1),inltx(1),pars,npars,&
                    0.0D+00,maxDose*1.1D+00,minvalue)
  ! estimate standard error of ED
  if(method==1) then
    ! 1) method 1, simple transformation
    averageErr=value/real(ndose,kind=8)
    lowltx=inltx(1)-sqrt((inltx(2))**2+averageErr)
    upltx =inltx(1)+sqrt((inltx(2))**2+averageErr)
    ! calculate low bounded ED, store it in lowup(1)
    call interpolate2(lowup(1),lowltx,pars,npars,&
                      0.0D+00,maxDose,minvalue)
    ! calculate up bounded ED, store it in lowup(2)
    call interpolate2(lowup(2),upltx,pars,npars,&
                      0.0D+00,maxDose*1.3D+00,minvalue)
    ! calculate standard error of ED with simple transformation
    outDose(2)=(lowup(2)-lowup(1))/2.0D+00
  else if(method==2) then
    ! 2) method 2, Monte Carlo iteration
    ! random seed
    seed=332571951
    ! initializing McCount numbers, Sum of Mc Dose, 
    ! Sum of the square of Mc Dose, i, they will be reused
    mccount=0
    sumdose=0.0D+00
    sumdose2=0.0D+00
    ! count number
    i=0
    ! Moto Carlo simulations
    do 
      ! set Ranpars to be pars
      Ranpars=pars
      ! generate random ltx values with mean=ltx(1,j), 
      ! sd=ltx(2,j), they will be store in ranltx
      do j=1, ndose
        call r8vec_normal(1,ltx(j,1),ltx(j,2),seed,ranltx(j))
      end do
      ! fitting x(Dose) .VS. random ltx
      if(npars==1) then
        call linearfit(Dose,ranltx,ndose,Ranpars,npars,.false.,&
                       Ranparserrors,Ranpredtval,Ranvalue,Ranerrorflag)
      else if(npars==2 .or. npars==3)  then
        call lmfit2(Dose,ranltx,ndose,Ranpars,npars,.false.,&
                    Ranparserrors,Ranpredtval,Ranvalue,&
                    lmtol,Ranerrorflag)
      end if
      ! error check
      if(Ranerrorflag(1)==123) then
        ! count really performed monte carlo numbers
        i=i+1
        ! generate random values signal value
        ! mcED(i) with mean=inltx(1), sd=inltx(2)
        call r8vec_normal(1,inltx(1),inltx(2),seed,mcED(i))
        ! interpolation using Ranpars and random values mcED(i), 
        ! calculated ED will store in the ith mcED
        call interpolate2(mcDose,mcED(i),Ranpars,npars,&
                          0.0D+00,maxDose*1.5D+00,minvalue)
        ! store the ith mcED
        mcED(i)=mcDose
        ! updating accumulated mcCount, SumDose, Sum of squared ED
        mccount=mccount+1
        sumdose=sumdose+mcDose
        sumdose2=sumdose2+mcDose**2
      end if
      ! terminate or not?
      if(i==motoiter) exit
    end do
    ! calculate standard error for ED
    outDose(2)=sqrt((real(mccount,kind=8)*sumdose2-sumdose**2)/&
                     real(mccount,kind=8)/real(mccount-1,kind=8))
  end if
  ! now return
  return
end subroutine calED2
