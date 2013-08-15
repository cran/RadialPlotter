subroutine aperr(ED,Error,nED,pars,sigma,&
                 spars,npars,tol,message)
!---------------------------------------------------------------------------------------------------
! subroutine aperr is used to approximate parameters' Std.Err
! for a minimum age model
! nED,              input:: integer, length of ED values (Error)
! npars,            input:: integer, number of parameters, either 3 or 4
! ED(nED),          input:: real values, ED values (unlogged)
! Error(nED),       input:: real values, absolute errors for ED values
! pars(npars),      input:: real values, parameters used for Std.Err calculation
! spars(npars),    output:: real values, estimated Std.Err
! tol,              input:: real values, the allowed maximum tolerance for the hessian to be non-sigular
! message,         output:: integer, error message generated during the calculation
! 
!     Author:: Peng Jun, 2013.07.24
!
! Dependence:: subroutine fdhessian; subroutine inverse; subroutine diag
!--------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::nED
  integer(kind=4),intent(in)::npars
  real   (kind=8),intent(in),dimension(nED)::ED
  real   (kind=8),intent(in),dimension(nED)::Error
  real   (kind=8),intent(in),dimension(npars)::pars
  real   (kind=8),intent(out),dimension(npars)::spars
  integer(kind=4),intent(out)::message
  real   (kind=8),intent(in)::sigma
  real   (kind=8),intent(in)::tol
  ! local variables
  real   (kind=8), parameter:: minAbsPar=0.0D+00
  real   (kind=8),dimension(npars,npars)::hessian
  real   (kind=8),dimension(npars)::grad
  integer(kind=4),dimension(4)::errorflag
  real   (kind=8),dimension(nED)::sED,sError
  real   (kind=8)::value
  integer(kind=4)::ifault
  ! return -99.0 if error appears
  spars=-99.0D+00
  message=0
  ! transform ED data to log-scale and add spread sigma value
  sError=sqrt((Error/ED)**2+sigma**2)
  sED=log(ED)
  ! call fdhessian to approximate hessian matrix using
  ! provided values of parameters
  call fdhessian(pars,sED,sError,npars,nED,tol, &
                 minAbsPar,hessian,grad,value,errorflag)
  ! error checking of fdhessian
  if(any(errorflag/=0))   then
    message=1
    return
  end if
  ! now inverse the hessian matrix to estimate the
  ! standard errors for calculated parameters 
  call inverse(hessian,npars,ifault,tol)
  ! check if error presents during inversing
  if(ifault/=0)  then
    message=1
    return
  end if
  ! extract diagnal elements of the 
  ! inversed hessian matrix
  call diag(hessian,npars,spars)
  ! check if any diagnal elment of inversed
  ! hessian is below zero
  if(any(spars<0.0D+00)) then
    message=1
    return
  end if
  spars=sqrt(spars)
  ! now return
  return
end subroutine aperr
