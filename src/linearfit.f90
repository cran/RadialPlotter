subroutine linearfit(xdat,ydat,ndat,pars,npars,calErr,&
                     stderror,predtval,value,info)
!--------------------------------------------------------------------------------------------------------------------
! linearfit is used to fit a dose-response cuerve in OSL dating of type linear  
!
! xdat(ndat),           input:: real values, the OSL dose date (De[Gy])
! ydat(ndat),           input:: real values, the standardlized OSL signal data (Lx/Tx)
! ndat,                 input:: integer, the length of dose or signal data
! pars(npars),         output:: real values, the initial guess pars, will be overwritten for output
! npars,                input:: integer, the length of pars
! calErr,               input:: logical, calculate pars' std.error or not
! stderror(npars),     output:: real values, estimated std.error for pars
! predtval(ndat),      output:: real values, fitted values correspond to ydat
! value,               output:: real value, sum of square of residuals
! info(2),             output:: integer, error message generated during the calling:
!                                1) if sum(xdat**2)=0 or sum(xxdat**2)=0, info(1)=123, else info(1)=0
!                                1) if pars' std.errors can be estimated, info(2)=0
!                                2) if error when calling lmhess to approximate model's hessian matrix, info(2)=1
!                                3) if hessian can not be inversed, info(2)=2
!                                4) if parameters' std.errors can not be estimated, info(2)=3
!
! Note:: if npars=1, a linear model of the form: y=ax will be fitted, where ndat>=1
!        if npars=2, a linear model of the form: y=ax+b will be fitted, where ndat>=2
!
! Author:: Peng Jun, 2013.09.21
!
! Dependence:: subroutine lmhess; subroutine inverse; subroutine diag
!------------------------------------------------------------------------------------------------------------------------
  
  implicit none
  integer(kind=4),intent(in)::ndat
  integer(kind=4),intent(in)::npars
  real   (kind=8),dimension(ndat),intent(in)::xdat,ydat
  real   (kind=8),dimension(npars),intent(out)::pars
  real   (kind=8),dimension(npars),intent(out)::stderror
  real   (kind=8),dimension(ndat),intent(out)::predtval
  real   (kind=8),intent(out)::value
  integer(kind=4),intent(out)::info(2)
  logical, intent(in):: calErr
  ! variables for subroutine lmhess
  real   (kind=8),parameter::lmtol=1.0D-07        ! maximum tolerance for singular matrix diagnosation 
  real   (kind=8),parameter::minAbsPar=0.0D+00    ! used in lmhess 
  real   (kind=8),dimension(npars,npars)::hessian ! hessian matrix obtained with finite-difference approximation
  real   (kind=8),dimension(npars)::gradient      ! gradient obtained with finite-difference approximation
  integer(kind=4),dimension(5)::hesserror         ! error message generated in subroutine lmhess
  real   (kind=8)::hessvalue                      ! value in subroutine lmhess
  ! variables for subroutine inverse
  integer(kind=4)::inverror                       ! error message generated in subroutine inverse
  ! local variables
  real   (kind=8),dimension(ndat)::xxdat
  integer(kind=4)::lmhessmodel
  !
  pars=-99.0D+00
  stderror=-99.0D+00
  predtval=-99.0D+00
  value=-99.0D+00
  info(1)=123
  info(2)=0
  ! calculate pars
  if(npars==1) then
    if(sum(xdat**2)<=lmtol) then
      info(1)=0
      return
    end if
    pars(1)=sum(xdat*ydat)/sum(xdat**2)
    predtval=pars(1)*xdat
    value=sum((ydat-predtval)**2)
  else if(npars==2) then
    if(sum(xdat**2)<=lmtol) then
      info(1)=0
      return
    end if
    xxdat=xdat-sum(xdat)/real(ndat,kind=8)
    pars(1)=sum(xxdat*ydat)/sum(xxdat**2)
    pars(2)=(sum(ydat)-sum(xdat)*pars(1))/real(ndat,kind=8)
    predtval=pars(1)*xdat+pars(2)
    value=sum((ydat-predtval)**2)
  end if
  if(calErr .eqv. .FALSE.)  return
  ! set lmhessmodel
  if(npars==1) then
    lmhessmodel=6
  else if(npars==2) then
    lmhessmodel=2
  end if
  ! calculate stderror
  call lmhess(pars,xdat,ydat,npars,ndat,lmtol,minAbsPar,&
              hessian,gradient,hessvalue,hesserror,lmhessmodel)
  ! check if any error appears when calling lmhess, 
  ! if so, reset info to be 1 and return
  if(any(hesserror/=0))  then
    info(2)=1
    return
  end if
  !
  ! set hessian to be its inverse 
  call inverse(hessian,npars,inverror,lmtol)
  ! check if any error appears when calling inverse, 
  ! if so, reset info to be 2 and return
  if(inverror==1)  then
    info(2)=2
    return
  end if
  ! extract diagnal elements from inversed hessian 
  ! matrix and storing it in arrary stderror
  call diag(hessian,npars,stderror)
  ! check if any elements in stderror is 
  ! below zero, if so, set info to be 3
  if(any(stderror<0)) info(2)=3
  ! rest stderror
  stderror=sqrt(stderror)
  ! terminate and return
  return
end subroutine linearfit
