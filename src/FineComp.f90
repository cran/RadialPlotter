subroutine FineComp(ED,Error,n,goodcomp,spreadsigma,&
                    maxiter,eps,maxcomp,BILI)
!-------------------------------------------------------------------------------------------------------
! FineComp is used to find out the best fitted component number of finite mixture age model,  
! it call subroutine FineED with various component numbers that range from 2 to maxcomp, and work
! out with the component number that given the minimum BIC value, note that the Equivalent
! Dose data must be unlogged
! 
! ED(n),input                 :: real values, the Equivalent Dose
! Error(n), input             :: real values, the assocaited absolute error of 
!                                Equivalent Dose                              
! n, input                    :: integer, the size of ED (or Error)
! goodcomp, output            :: integer, the estimated finest component number
! spreadsigma, input          :: real value, the spread dispersion to the relative 
!                                error of Equivalent Dose
! maxiter, input              :: integer, the maximum iterative number allowed
! eps, input                  :: real value, the maximum tolerance for converge
! maxcomp, input              :: the maximum allowed component number for loopping
! BILI(maxcomp-1,2), output   :: the BIC and maxlik values for component number ranging
!                                from 2 to maxcomp
!
! Author:: Peng Jun, 2013.03.03, revised in 2013.03.04
!
! References :: Galbraith RF, 1988. Graphical Display of Estimates Having Differing 
!               Standard Errors.Techno-metrics, 30, page 271-281.
!               Galbraith, RF, Green PF, 1990. Estimating the component ages in a 
!               finite mixture. Nuclear Tracks and Radiation Measurements, 17, page 197-206.
!               
! Dependence  :: subroutine FineED
!
!------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::n
  integer(kind=4),intent(out)::goodcomp
  integer(kind=4),intent(in)::maxiter
  integer(kind=4),intent(in)::maxcomp
  real   (kind=8),dimension(n),intent(in)::ED
  real   (kind=8),dimension(n),intent(in)::Error
  real   (kind=8),intent(in)::spreadsigma 
  real   (kind=8),intent(in)::eps 
  real   (kind=8),dimension(maxcomp-1,2),intent(out)::BILI
  ! local variables
  real(kind=8)::minBIC
  real(kind=8)::maxlik
  real(kind=8)::BIC
  real(kind=8),allocatable::pars(:,:) 
  integer(kind=4)::i,j
  !
  minBIC=1.0D+30
  goodcomp=0
  !
  do i=2,maxcomp
    !
    allocate(pars(1:2,1:i))
    !
    call FineED(ED,Error,n,i,spreadsigma,pars,&
                maxlik,BIC,maxiter,eps)
    BILI(i-1,1)=BIC
    BILI(i-1,2)=maxlik
    !
    if(BIC<minBIC)  then
      minBIC=BIC
      goodcomp=i
    end if
    deallocate(pars)
  end do
  !
  return
  !
end subroutine FineComp
