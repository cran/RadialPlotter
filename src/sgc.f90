subroutine sgc(dose,ltx,ndose,inltx,outDose,nout,pars,predtval,&
               parserrors,npars,value,errorflag,method,motoiter)
!-----------------------------------------------------------------------------------------------------------------------------
! sgc is a subroutine used for analyzing Standardlised Grow Curve method in OSL dating.
!
! ndose,               input:: integer, length of Dose(Redose1, Redose2,...) or Ltx(Lx1/Tx1, Lx2/Tx2,...)
! nout,                input:: integer, length of standardlized OSL signal that to be analyzed
! npars,               input:: integer, model used for fitting the dose-response curve:
!                              if npars=2, a linear model of the form: y=ax+b will be fitted, where ndat>=2
!                              if npars=3, a Exponential model of the form: y=a(1-exp(-bx))+c  will be fitted, where ndat>=3
!                              if npars=4, a linear+Exponential model of the form: y=a(1-exp(-bx))+cx+d will be fitted, where ndat>=4
! dose(ndose),         input:: real values, Redose data for dose-response curve
! ltx(ndose),          input:: real values, standardlized OSL signal data for dose-response curve
! inltx(nout),         input:: real values, standardlized OSL signal values that used for calculating dose values
! outDose(nout,2),    output:: real values, calculated dose and standard errors that correlated to standardlized OSL signal
! pars(npars),  input/output:: real values, initial guess characteristical parameters for dose-response curve, being overwritten as ouput
! predtval(ndose),    output:: real values, fitted values that correspond to ltx
! parserrors(npars),  output:: real values, parameters' standard errors
! value,              output:: real value, sum of the square of residual
! motoiter,           input :: integer, numbers of Moto Carlo simulations when estimating standard error for each dose value
! method,              input:: integer, method for calculating dose's standard error, 1 for simple method; 2 for Moto Carlo method
! errorflag(2),       output:: integer, error message generated when calling subroutine SGC:
!                              1) error appears when calling subroutine lmder1, errorflag(1) will be one in (0,4,5,6,7,8,9,10), 
!                                 else errorflag(1)=123
!                              2) error appears when attempt to calculate parameters' standard errors, errorflag(2)=1, else
!                                 errorflag(2) will be 0
!   interpolate(Dose,ltx,pars,npars,lowb,upb,typ,value)
!
! Author:: Peng Jun, 2013.05.26, revised in 2013.05.29
!
! Dependence:: subroutine lmfit1; subroutine interpolate; subroutine r8vec_normal
! 
! References:: Roberts,H.M. and Duller,G.A.T., 2004. Standardised growth curves for optical dating of sediment
!              using multiple-grain aliquots. Radiation Measurements 38, pp. 241-252.
!
!              Duller, G.A.T., 2007. Assessing the error on equivalent dose estimates derived from single aliquot
!              regenerative dose measurements. Ancient TL 25, pp. 15-24.
!---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::ndose
  integer(kind=4),intent(in)::nout
  integer(kind=4),intent(in)::npars
  integer(kind=4),intent(in)::motoiter
  integer(kind=4),intent(in)::method
  real   (kind=8),intent(out)::value
  integer(kind=4),intent(out)::errorflag(2)
  real   (kind=8),dimension(ndose),intent(in)::dose,ltx
  real   (kind=8),dimension(nout,2),intent(in)::inltx
  real   (kind=8),dimension(nout,2),intent(out)::outDose
  real   (kind=8),dimension(npars),intent(inout)::pars
  real   (kind=8),dimension(npars),intent(out)::parserrors
  real   (kind=8),dimension(ndose),intent(out)::predtval
  !
  ! variables for subroutine lmfit1
  real(kind=8),parameter::tol=1.0D-07       ! toleracne for stopping iterations in subroutine lmfit1
  !
  ! variables for subroutine interpolate
  real   (kind=8)::mcDose                   ! Dose that calcualted using Moto Carlo method
  real   (kind=8)::lowb, upb                ! interval from which the interpolation to take place
  real   (kind=8)::minvalue,cminvalue       ! the minimized value by interpolating
  real   (kind=8)::charDose                 ! a Dose value that given minimized value
  logical:: typ                             ! type of interpolation
  real   (kind=8)::lowltx,upltx             ! low and up boundary ltx value
  ! local variables
  integer(kind=4)::i,j,mccount              ! looping counter
  real(kind=8),dimension(nout,2)::lowup     ! low and up boundary dose values
  real(kind=8),dimension(motoiter)::moto    ! random standardlized OSL signal generated using Moto Carlo method
  real(kind=8)::sumdose,sumdose2            ! sum of dose and sum of the square of dose for each Moto Carlo based dose value
  integer(kind=4)::seed                     ! seed for normal distribution random number generation
  real(kind=8)::averageErr                  ! average value for the sum of the squared residual errors
  !
  !****************************111*******************************************
  ! call lmfit1 to get the characterized parameters for dose-response 
  ! curve, pars will be reset, also return parserrors,predtval,value, errorflag
   call lmfit1(dose,ltx,ndose,pars,npars,&
               parserrors,predtval,value,tol,errorflag)
  !
  !*****************************222******************************************
  ! using calculated characterized parameters to do interpoltating to
  ! calculate Dose values from the Dose-Response Curve
  !
  ! fisrtly, calculate the maximum ltx value using the characterized 
  ! dose-response curve. 
  ! set lowb and upb 
  lowb=0.0D+00                      ! lowb is set to be zero
  upb=maxval(Dose)*2.0D+00          ! upb is set to be mutiple of the maximum value in dose
  call interpolate(charDose,0.0D+00,pars,npars,lowb,upb,.false.,minvalue)
  ! then, calculate Dose value using corresponded ltx value 
  do i=1, nout
    ! check if a Dose values can be calculated
    if(inltx(i,1)<= -minvalue) then
      call interpolate(outDose(i,1),inltx(i,1),pars,npars,lowb,upb,.true.,cminvalue)
    else 
      ! if can not be estimated, set it to be -99.0
      outDose(i,1)=-99.0D+00
    end if
    !
  end do
  !
  !*****************************333********************************************
  ! calculate standard errors for equivalent dose from the Dose-Response Curve
  !
  ! if standard errors for Lx/Tx can not be provided,
  ! Equivalent Doses' standard errors can not be calculated,
  ! just return equivalent dose values
  if(all(inltx(:,2)<=1.0D-13) )  then
    outDose(:,2)=-99.0D+00
    return
  else
    ! if standard errors for Lx/Tx can be obtained, then
    ! calculating Equivalent Doses' standard errors.
    ! method=1 will using simple methods to estimated
    ! Equivalent Doses' standard errors
    if(method==1) then
      ! calculate average value for the sum of the squared residual errors,
      ! which will be bundled to the ith standard errors of Lx/Tx
      averageErr=value/real(ndose,kind=8)
      do i=1,nout
        ! check if low boundary Dose values can be calculated
        lowltx=inltx(i,1)-sqrt((inltx(i,2))**2+averageErr)
        if(lowltx<= -minvalue) then
          ! calculate low boundary Dose values
          call interpolate(lowup(i,1),lowltx,pars,npars,lowb,upb,.true.,cminvalue)
        else 
          ! if can not be estimated, set it to be -99.0
          lowup(i,1)=-99.0D+00
        end if      
        ! check if up boundary Dose values can be calculated
        upltx=inltx(i,1)+sqrt((inltx(i,2))**2+averageErr)
        if(upltx<= -minvalue) then
          ! calculate up boundary Dose values
          call interpolate(lowup(i,2),upltx,pars,npars,lowb,upb,.true.,cminvalue)
        else 
          ! if can not be estimated, set it to be -99.0
          lowup(i,2)=-99.0D+00
        end if 
      end do
      ! estimate Doses's standard errors using simple transformation
      outDose(:,2)=(lowup(:,2)-lowup(:,1))/2.0D+00
    !
    ! if method=2, Monte Carlo method will be used in standard error estimation
    else if(method==2) then 
      ! seed for normal distribution number generation
      seed=532571951
      ! for i in 1 to maximum number of inltx do
      do i=1,nout
        ! generate normal distribution random number
        ! with mean=inltx(i,1) and stdev=inltx(i,2)
        call r8vec_normal (motoiter, inltx(i,1), inltx(i,2), seed, moto)
        ! set sum of dose value and sum of the square of dose value, also a counter
        mccount=0
        sumdose=0.0D+00
        sumdose2=0.0D+00
        ! add up sumdose and sumdose2 to calculate the standard error for each dose value
        do j=1,motoiter
          if(moto(j)<=-minvalue) then         
            call interpolate(mcDose,moto(j),pars,npars,lowb,upb,.true.,cminvalue)  
            mccount=mccount+1
            sumdose=sumdose+mcDose
            sumdose2=sumdose2+mcDose**2
          end if
        end do
        ! calculate standard error for the ith dose value, store it 
        ! in the ith row of the second column of matrix outDose
        outDose(i,2)=sqrt( (real(mccount,kind=8)*sumdose2-sumdose**2)/ &
                     real(mccount,kind=8)/real(mccount-1,kind=8) )
      end do
    end if
  end if
  ! at this point all values needed for output
  ! have been calculated, just return
  return
end subroutine sgc
