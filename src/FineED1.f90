subroutine FineED1(ED,Error,n,ncomp,spreadsigma,pars,&
                   spars,maxlik,BIC,maxiter,eps,tol,errorflag)
!-------------------------------------------------------------------------------------------------------
! FineED is used to perfect the estimated finite mixture age model based Equivalent Dose,
! it first change the data to log-scale and initialize the parameters automatically with various
! ranges, then call subroutine FMMED to work out the parameter that gives the minimum BIC value
! it different from subroutine FineED that it also calculates the parameters' standard error, it also
! allows the perform of Central Age Model analysis by specifying ncomp=1
!
! ED(n),input                 :: real values, the Equivalent Dose, must be unlogged
! Error(n), input             :: real values, the assocaited absolute error of 
!                                Equivalent Dose                          
! n, input                    :: integer, the size of ED (or Error)
! ncomp, input                :: integer, the component number
! spreadsigma, input          :: real value, the spread dispersion to the relative 
!                                error of logged Equivalent Dose
! pars(2,ncomp), output       :: real values, the estimated parameters
! spars(2,ncomp), output      :: real values, the standard error for estimated parameters
! maxlik, output              :: real value, the maximum logged-likelihood value
! BIC, output                 :: real value, the BIC value
! maxiter, input              :: integer, the maximum iterative number allowed
! eps, input                  :: real value, the maximum tolerance for stopping the
!                                iterative process
!
! Author:: Peng Jun, 2013.03.05
!
! References :: Galbraith RF, 1988. Graphical Display of Estimates Having Differing 
!               Standard Errors.Techno-metrics, 30, page 271-281.
!
!               Galbraith, RF, Green PF, 1990. Estimating the component ages in a 
!               finite mixture. Nuclear Tracks and Radiation Measurements, 17, page 197-206.
!               
! Dependence  :: subroutine FMMED; subroutine kmeans; 
!                subroutine CAM; subroutine AppCovar
!
!------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::n
  integer(kind=4),intent(in)::ncomp
  integer(kind=4),intent(in)::maxiter
  integer(kind=4),intent(out)::errorflag
  real   (kind=8),dimension(n),intent(in)::ED
  real   (kind=8),dimension(n),intent(in)::Error
  real   (kind=8),dimension(2,ncomp),intent(out)::pars
  real   (kind=8),dimension(2,ncomp),intent(out)::spars
  real   (kind=8),intent(in)::spreadsigma 
  real   (kind=8),intent(out)::maxlik
  real   (kind=8),intent(out)::BIC
  real   (kind=8),intent(in)::eps 
  real   (kind=8),intent(in)::tol
  !
  !local variables for subroutine kmeans
  integer(kind=4),parameter::iter=20
  integer(kind=4),parameter::nstart=100
  integer(kind=4)::belo(n),clusp(ncomp)
  real   (kind=8)::energ(ncomp),clust(ncomp)
  integer(kind=4)::ifault
  !
  ! local variables
  integer(kind=4)::i,j
  real(kind=8)::minBIC
  real(kind=8)::sMaxlik
  real(kind=8),dimension(2,ncomp)::cpars,kpars
  real(kind=8),dimension(n)::sED,sError
  !
  sError=Error/ED
  sED=log(ED)
  !
  if (ncomp==1)  then 
    errorflag=0
    call CAM(ED,Error,n,spreadsigma,maxiter,&
             eps,pars,spars,maxlik,BIC)
    return
  else 
    ! initialize cpars
    cpars(1,:)=pars(1,:)
    do j=1,ncomp
      cpars(2,j)=minval(sED)+(maxval(sED)-minval(sED))*real(j-1)/real(ncomp-1)
    end do
    !
    minBIC=1.0D+30
    !
    do i=1,3
      !
      pars(1,:)=1.0D+00/real(ncomp)
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
    call kmeans(1,n,ncomp,iter,nstart, &
                sED,belo,clust,clusp,energ,ifault)
    !
    kpars(1,:)=1.0D+00/real(ncomp)
    kpars(2,:)=clust
    !
    call FMMED(sED,sError,n,ncomp,spreadsigma,kpars,&
               maxlik,BIC,maxiter,eps)
    !
    if(BIC<minBIC)  then
      cpars=kpars
      minBIC=BIC
      sMaxlik=maxlik
    end if
    !
    pars=cpars
    BIC=minBIC
    maxlik=sMaxlik
    !
    ! approximate standard errors
    call AppCovar(pars,ncomp,ED,Error,n, &
                  spreadsigma,spars,errorflag,tol)
    !
    ! change equivalent dose to un-logged scale
    pars(2,:)=exp(pars(2,:))
    spars(2,:)=spars(2,:)*pars(2,:)
    !
    return
  end if
  !
end subroutine FineED1
