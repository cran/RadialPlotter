subroutine CAM(ED,Error,n,spreadsigma,maxiter,&
               eps,pars,spars,maxlik,BIC)
!------------------------------------------------------------------------------------------
! CAM is used to do Central Age Model based Equivalent Dose calculation,
! this version only used logged equivalent Dose for analyzing
!
! ED(n),input                   :: real values, the logged Equivalent Dose 
! Error(n), input               :: real values, the assocaited relative error 
!                                  for Equivalent Dose
! n, input                      :: integer, the size of ED (or Error)
! spreadsigma, input,           :: real value, the spread dispersion to the relative 
!                                  error of logged Equivalent Dose
! maxiter, input                :: integer, the maximum iterative number allowed
! eps, input                    :: real value, the maximum tolerance for stop the loop
! pars(2,1), output             :: real values, the estimated parameters
! spars(2,1), output            :: real values, the standard errors for estimated parameters
! maxlik, output                :: real value, the maximum logged-likelihood value
! BIC,    output                :: real value, the BIC value
!
! Author:: Peng Jun, 2013.03.04, revised in 2013.03.05
!
! Dependence:: No
!
! References:: Galbraith RF, 1988. Graphical Display of Estimates Having Differing 
!              Standard Errors.Techno-metrics, 30, page 271-281.
!              Galbraith, RF, Green PF, 1990. Estimating the component ages in a 
!              finite mixture. Nuclear Tracks and Radiation Measurements, 17, page 197-206.
! 
!--------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::n
  integer(kind=4),intent(in)::maxiter
  real   (kind=8),intent(in)::eps
  real   (kind=8),intent(in)::spreadsigma
  real   (kind=8),intent(out)::maxlik
  real   (kind=8),intent(out)::BIC
  real   (kind=8),dimension(n),intent(in)::ED
  real   (kind=8),dimension(n),intent(in)::Error
  real   (kind=8),dimension(2,1),intent(out)::pars
  real   (kind=8),dimension(2,1),intent(out)::spars
  ! local variables
  real   (kind=8),dimension(n)::z,sz,wz
  real   (kind=8),dimension(2,1)::newpars
  integer(kind=4)::i
  ! change data to log-scale
  z=log(ED)
  sz=sqrt((Error/ED)**2+spreadsigma**2)
  ! guess initial sigma
  pars(1,1)=0.1D+00
  !
  wz=1.0D+00/((pars(1,1))**2+sz**2)
  !
  pars(2,1)=sum(z*wz)/sum(wz)
  !
  ! do the major loop
  do i=1,maxiter
    !
    newpars(2,1)=sum(z*wz)/sum(wz)
    newpars(1,1)=pars(1,1)*sqrt(sum((wz**2)*(z-newpars(2,1))**2/sum(wz)))
    wz=1.0D+00/((newpars(1,1))**2+sz**2)
    !
    if (abs(pars(1,1)-newpars(1,1))+&
        abs(pars(2,1)-newpars(2,1))<eps)  then
       goto 100
     else
       pars=newpars
     end if
     !
   end do
  !
  100 maxlik=0.5D+00*sum(log(wz))-0.5D+00*sum(wz*(z-pars(2,1))**2)
  BIC=-2.0D+00*maxlik-2.0D+00*log(real(n))
  pars(2,1)=exp(pars(2,1))
  spars(2,1)=pars(2,1)/sqrt(sum(wz))
  spars(1,1)=1.0D+00/sqrt(2.0D+00*pars(1,1)*sum(wz**2))
  !
  return
  !
end subroutine CAM  
