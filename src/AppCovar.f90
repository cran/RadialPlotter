subroutine AppCovar(pars,ncomp,ED,Error,n,&
                    spreadsigma,ParError,iferror,tol)
!------------------------------------------------------------------------------------------------------
! AppCovar is used to estimate the standard error of parameters using
! covariance matrix approximation, note that the Equivalent Dose data
! must be of unlogged
! 
! pars(2,ncomp),      input:: real values, the parameters estimated by subroutine FMMED
! ncomp,              input:: integer, the component number
! ED(n),              input:: real values, the Equivalent Dose
! Error(n),           input:: real values, the absolute error for Equivalent Dose
! n      ,            input:: integer, the size of ED (or Error)
! spreadsigma,        input:: real value, the spread dispersion to the error of Equivalent Dose
! ParError(2,ncomp), output:: the estimated standard error by covariance matrix approximation
! iferror,           output:: integer, the error message when inversing the covariance matrix, 
!                             0 for a matrix that is both non-singular and have non-minus  
!                             elements in its iverse matrix, otherwise 1
! tol,                input:: real value, diagnosing A to be a singular matrix if max(abs(A)) 
!                             is smaller than tol
!
! Author:: Peng Jun, 2013.02.27
!
! References :: Galbraith RF, 1988. Graphical Display of Estimates Having Differing Standard Errors. 
!               Techno-metrics, 30, page 271-281.
!
! Dependence:: subroutine inverse, subroutine diag
!
!----------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::ncomp
  integer(kind=4),intent(in)::n
  integer(kind=4),intent(out)::iferror
  real   (kind=8),intent(in)::tol
  real   (kind=8),intent(in)::spreadsigma
  real   (kind=8),dimension(2,ncomp),intent(in)::pars
  real   (kind=8),dimension(n),intent(in)::ED
  real   (kind=8),dimension(n),intent(in)::Error
  real   (kind=8),dimension(2,ncomp),intent(out)::ParError
  ! local variables
  real(kind=8),dimension(n)::z
  real(kind=8),dimension(n)::s
  real(kind=8),dimension(n,ncomp)::p
  real(kind=8),dimension(n,ncomp)::pf
  real(kind=8),dimension(n)::rspf
  real(kind=8),dimension(ncomp,ncomp)::II
  real(kind=8),dimension(ncomp-1,ncomp-1)::A
  real(kind=8),dimension(ncomp-1,ncomp)::B
  real(kind=8),dimension(ncomp,ncomp)::C
  real(kind=8),dimension(n)::w
  real(kind=8),dimension(n,ncomp)::AA
  real(kind=8),dimension(n,ncomp)::BB
  real(kind=8),dimension(2*ncomp-1)::diagcovar
  real(kind=8),dimension(2*ncomp-1,2*ncomp-1)::covar
  integer(kind=4)::i,j
  !
  z=log(ED)
  s=Error/ED
  w=1.0D+00/(spreadsigma**2+s**2)
  !
  AA=0.0D+00
  BB=0.0D+00
  A=0.0D+00
  B=0.0D+00
  C=0.0D+00
  !
  do i=1,ncomp
    pf(:,i)=sqrt(w)*exp(-0.5D+00*w*(z-pars(2,i))**2)*pars(1,i)
  end do
  !
  rspf=sum(pf,2)
  !
  do j=1, n
    p(j,:)=pf(j,:)/rspf(j)
  end do
  !
  II=0.0D+00
  !
  do i=1,ncomp
    AA(:,i)=w*(z-pars(2,i))
    BB(:,i)=-w+(w*(z-pars(2,i)))**2
    II(i,i)=1.0D+00
  end do
  !
  do i=1,ncomp-1
    do j=1,ncomp-1
      A(i,j)=sum( (p(:,i)/pars(1,i)-p(:,ncomp)/pars(1,ncomp))* &
                  (p(:,j)/pars(1,j)-p(:,ncomp)/pars(1,ncomp))  )
    end do
  end do
  !
  do i=1,ncomp-1
    do j=1,ncomp
      B(i,j)=sum( p(:,j)*AA(:,j)*(p(:,i)/pars(1,i)- &
                  p(:,ncomp)/pars(1,ncomp)- &
                  II(i,j)/pars(1,i)+ &
                  II(ncomp,j)/pars(1,ncomp)) )
    end do
  end do
  !
  do i=1,ncomp
    do j=1,ncomp
      C(i,j)=sum( p(:,i)*p(:,j)*AA(:,i)*AA(:,j)- &
                  II(i,j)*BB(:,i)*p(:,i) )
    end do
  end do
  !
  covar(1:ncomp-1,1:ncomp-1)=A
  covar(1:ncomp-1,ncomp:)=B
  covar(ncomp:,1:ncomp-1)=transpose(B)
  covar(ncomp:,ncomp:)=C
  !
  call inverse(covar,2*ncomp-1,iferror,tol)
  !
  call diag(covar,2*ncomp-1,diagcovar)
  !
  if(any(diagcovar<0.0D+00)) iferror=1
  !
  ParError(1,1:ncomp-1)=sqrt(diagcovar(1:ncomp-1))
  ParError(1,ncomp)=sqrt(sum(covar(1:ncomp-1,1:ncomp-1)))
  !
  ParError(2,:)=sqrt(diagcovar(ncomp:))
  !
  return
  !
end subroutine AppCovar     
