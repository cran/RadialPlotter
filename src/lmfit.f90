subroutine lmfit(xdat,ydat,ndat,pars,npars,stderror,predtval,value,tol,info)
!--------------------------------------------------------------------------------------------------------------------------------
! specifying a fitting model, lmfit will fit the model using Levenberg-Marquadt method
! model=1 has the formula I(t)=a1*exp(-b1*t)+a2*exp(-b2*t)+...+ak*exp(-bk*t), K=1:7
!
! ndat,                input:: integer, length of xdat or y dat.
! npars,               input:: integer, dimension of parameters (or length of pars).
! xdat(ndat),          input:: real values, xdat (or independent variables x).
! ydat(ndat),          input:: real values, ydat (or dependent variables y).
! pars(npars),  input/output:: real values, initial guess values for parameters to be optimized,
!                              overwritten to be final results for output.
! stderror(npars),    output:: real values, estimated standard errors for parameters
! predtval(ndat),     output:: real values, fited values correspond to ydat
! value,              output:: real value, final total residual error
! tol,                 input:: real value,  Termination occurs when the algorithm estimates either that the relative error in the  
!                              sum of squares is at most TOL or that the relative error between X and the solution is at most TOL.
! info,               output:: error message generated during the calculation:
!                              1) successful work, info=1
!                              2) error when calling lmder1, info=2
!                              3) at least one estimated parameter is below zero, info=3
!                              4) error when calling lmhess to approximate model's hessian matrix, info=4
!                              5) error when attempt to inverse the approximated hessian matrix, info=5
!                              6) error when attempt to calculate pars' standard errors, info=6
!
! Author:: Peng Jun, 2013.05.21, revised in 2013.05.22
!
! Dependence:: subroutine lmfunc; subroutine lmder1; subroutine lmhess; subroutine inverse; subroutine diag.
!---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::ndat
  integer(kind=4),intent(in)::npars
  real   (kind=8),dimension(ndat),intent(in)::xdat,ydat
  real   (kind=8),dimension(npars),intent(inout)::pars
  real   (kind=8),dimension(npars),intent(out)::stderror
  real   (kind=8),dimension(ndat),intent(out)::predtval
  real   (kind=8),intent(in)::tol
  real   (kind=8),intent(out)::value
  integer(kind=4),intent(out)::info
  !
  ! variables for subroutine lmhess
  real   (kind=8),parameter::lmtol=1.0D-07        ! used for singular matrix diagnose in subroutine lmhess and subroutine inverse
  real   (kind=8),parameter::minAbsPar=0.0D+00    ! for lmhess use
  real   (kind=8),dimension(npars,npars)::hessian ! hessian matrix by finite-difference approximation
  real   (kind=8),dimension(npars)::gradient      ! gradient by finite-difference approximation
  integer(kind=4),dimension(5)::hesserror         ! error message generated in lmhess
  ! variables for subroutine inverse
  integer(kind=4)::inverror                       ! for singular matrix diagnose in subroutine inverse
  ! local variables
  integer(kind=4)::ldfjac,iflag
  real   (kind=8),dimension(ndat)::fvec           ! fitted residual in vector form
  real   (kind=8),dimension(ndat,npars)::fjac     ! jacobian matrix
  !
  ldfjac=ndat
  !
  ! default stderrors if error appears
  stderror=-99.0D+00
  ! now specifying iflag to be 1 to calculate fvec
  iflag=1
  !
  ! using initial pars to caculate fvec that will be used in subroutine lmder1
  call lmfunc(xdat,ydat,ndat,npars,pars,fvec,fjac,ldfjac,iflag)
  !
  ! optimizing initial pars using Levenberg-Marquadt method
  ! and return pars and info for output
  call lmder1(lmfunc,ndat,npars,pars,fvec,fjac,ldfjac,tol,info,xdat,ydat)
  ! calculate sum of squre of residual 
  value=sum( (fvec)**2 )
  ! calculate fitted values correspond to ydat
  predtval=fvec+ydat
  !
  ! check if any error appears when calling lmder1, and reset info 
  ! to 2 and return if so, else reset info to 1 and continue
  if(info==1 .or. info==2 .or. info==3 ) then
    info=1
  else 
    info=2
    return
  end if  
  !
  ! check if any estimated parameter is below zero,
  ! if so, reset info to be 3 and return
  if(any(pars<0.0D+00) ) then
    info=3
    return
  end if
  ! estimate pars' standard errors
  call lmhess(pars,xdat,ydat,npars,ndat,lmtol,minAbsPar, &
              hessian,gradient,value,hesserror,1)
  ! reset value
  value=value**2
  ! check if any error appears when calling lmhess, 
  ! if so, reset info to be 4 and return
  if( any(hesserror/=0) )  then
    info=4
    return
  end if
  !
  ! set hessian to be its inverse 
  call inverse(hessian,npars,inverror,lmtol)
  ! check if any error appears when calling inverse, 
  ! if so, reset info to be 5 and return
  if(inverror==1)  then
    info=5
    return
  end if
  ! extract diagnal elements from inversed hessian 
  ! matrix and storing it in arrary stderror
  call diag(hessian,npars,stderror)
  ! check if any elements in stderror is 
  ! below zero, if so, set info to be 6
  if( any(stderror<0) ) info=6
  ! rest stderror
  stderror=sqrt(stderror)
  return
end subroutine lmfit