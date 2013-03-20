subroutine FineED(ED,Error,n,ncomp,spreadsigma,&
                  pars,maxlik,BIC,maxiter,eps)
!-------------------------------------------------------------------------------------------------------
! FineED is used to perfect the estimated finite mixture age model based Equivalent Dose,
! it first change the data to log-scale and initialize the parameters automatically with various
! ranges, then call subroutine FMMED to work out the parameter that gives the minimum BIC value
!
! ED(n),input                 :: real values, the Equivalent Dose, must be unlogged
! Error(n), input             :: real values, the assocaited absolute error of 
!                                Equivalent Dose                          
! n, input                    :: integer, the size of ED (or Error)
! ncomp, input                :: integer, the component number
! spreadsigma, input          :: real value, the spread dispersion to the relative 
!                                error of logged Equivalent Dose
! pars(2,ncomp), input        :: real values, the estimated parameters
! maxlik, output              :: real value, the maximum logged-likelihood value
! BIC, output                 :: real value, the BIC value
! maxiter, input              :: integer, the maximum iterative number allowed
! eps, input                  :: real value, the maximum tolerance for stopping the
!                                iterative process
!
! Author:: Peng Jun, 2013.03.03
!
! References :: Galbraith RF, 1988. Graphical Display of Estimates Having Differing 
!               Standard Errors.Techno-metrics, 30, page 271-281.
!               Galbraith, RF, Green PF, 1990. Estimating the component ages in a 
!               finite mixture. Nuclear Tracks and Radiation Measurements, 17, page 197-206.
!               
! Dependence  :: subroutine FMMED
!
!------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::n
  integer(kind=4),intent(in)::ncomp
  integer(kind=4),intent(in)::maxiter
  real   (kind=8),dimension(n),intent(in)::ED
  real   (kind=8),dimension(n),intent(in)::Error
  real   (kind=8),dimension(2,ncomp),intent(out)::pars
  real   (kind=8),intent(in)::spreadsigma 
  real   (kind=8),intent(out)::maxlik
  real   (kind=8),intent(out)::BIC
  real   (kind=8),intent(in)::eps 
  ! local variables
  integer(kind=4)::i,j
  real(kind=8)::minBIC
  real(kind=8)::sMaxlik
  real(kind=8),dimension(2,ncomp)::cpars
  real(kind=8),dimension(n)::sED,sError
  !
  sError=Error/ED
  sED=log(ED)
  minBIC=1.e30
  pars(1,:)=1.0/real(ncomp)
  !
  do i=1,3
    !
    do j=1,ncomp
      pars(2,j)=minval(sED)+(maxval(sED)-minval(sED))*real(j)/real(ncomp+i-2)
    end do
    !
    call FMMED(sED,sError,n,ncomp,spreadsigma,pars,&
               maxlik,BIC,maxiter,eps)
    !
    if(BIC<minBIC)  then
      cpars=pars
      minBIC=BIC
      sMaxlik=maxlik
    end if
    !
  end do
  !
  pars=cpars
  BIC=minBIC
  maxlik=sMaxlik
  !
  return
  !
end subroutine FineED
