subroutine FMMED(ED,Error,n,ncomp,spreadsigma,&
                 pars,maxlik,BIC,maxiter,eps)
!-------------------------------------------------------------------------------------------------------
! FMMED is used to do finite mixture age analysis by specify the initial guess parameters, then
! it converge to a local minimum value by iterative Maximum Likelihood Method, both logged and
! unlogged Equivalent Dose can be used for analyzing
!
! ED(n),input                   :: real values, the Equivalent Dose (or logged Equivalent Dose) 
! Error(n), input               :: real values, the assocaited absolute error (or relative error)                             
! n, input                      :: integer, the size of Equivalent Dose (or Error)
! ncomp, input                  :: integer, the component number
! spreadsigma, input            :: real value, the spread dispersion to the relative 
!                                  error of logged Equivalent Dose
! pars(2,ncomp), input/output   :: real values, the estimated parameters
! maxlik, output                :: real value, the maximum logged-likelihood value
! BIC, output                   :: real value, the BIC value
! maxiter, input                :: integer, the maximum iterative number allowed
! eps, input                    :: real value, the maximum tolerance for stop the iterative loop
!
! Author:: Peng Jun, 2013.02.26, revised in 2013.03.03
!
! References :: Galbraith RF, 1988. Graphical Display of Estimates Having Differing 
!               Standard Errors.Techno-metrics, 30, page 271-281.
!               Galbraith, RF, Green PF, 1990. Estimating the component ages in a 
!               finite mixture. Nuclear Tracks and Radiation Measurements, 17, page 197-206.
!               
! Dependence  :: No
!
!------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::n
  integer(kind=4),intent(in)::ncomp
  integer(kind=4),intent(in)::maxiter
  real   (kind=8),dimension(n),intent(in)::ED
  real   (kind=8),dimension(n),intent(in)::Error
  real   (kind=8),dimension(2,ncomp),intent(inout)::pars
  real   (kind=8),intent(in)::spreadsigma 
  real   (kind=8),intent(out)::maxlik
  real   (kind=8),intent(out)::BIC
  real   (kind=8),intent(in)::eps 
  !
  ! local variables
  real(kind=8),dimension(n)::rspf
  real(kind=8),dimension(n)::w
  real(kind=8),dimension(2,ncomp)::newpars
  real(kind=8),dimension(n,ncomp)::f
  real(kind=8),dimension(n,ncomp)::p,pp
  real(kind=8),dimension(n,ncomp)::wp
  real(kind=8),dimension(n,ncomp)::pf 
  real(kind=8),parameter::PI=3.141592653589793238462643383279502884197
  integer(kind=4)::i,j
  !
  ! set the weight value
  w=1.0/(spreadsigma**2+Error**2)
  !
  ! start the major loop
  do i=1,maxiter
    !
    do j=1,ncomp
      f(:,j)=sqrt(w)*exp(-0.5*w*(ED-pars(2,j))**2)
      pf(:,j)=pars(1,j)*f(:,j)
    end do
    !
    rspf=sum(pf,2)
    !
    do j=1,n
      pp(j,:)=pf(j,:)/rspf(j)
      wp(j,:)=w(j)*pp(j,:)
      p(j,:)=pp(j,:)*w(j)*ED(j)
    end do
    !
    ! updating characterized ED values
    ! and proportion values
    newpars(2,:)=sum(p,1)/sum(wp,1)
    newpars(1,:)=sum(pp,1)/real(n)
    !
    ! calculate the maximum logged-likelihood value
    ! and BIC value
    maxlik=sum(log(1.0/sqrt(2.0*PI)*sum(pf,2)))
    BIC=-2.0*maxlik+(2.0*real(ncomp)-1.0)*log(real(n))
    !
    ! test for converge
    if( sum(abs(pars(2,:)-newpars(2,:)))+&
        sum(abs(pars(1,:)-newpars(1,:)))<eps )  then
      return
    else
      pars(2,:)=newpars(2,:)
      pars(1,:)=newpars(1,:)
    end if
  !
  end do
  !
  return
  !
end subroutine FMMED
