subroutine fdhessian(pars,ED,Error,npars,nED,tol,minAbsPar, &
                     hessian,gradient,value,errorflag)
!--------------------------------------------------------------------------------------------------
! subroutine fdhessian is used to calculate the gradient, hessian matrix of a 
! given function contained in its inner, that is fun34 in MAM models, using 
! finite-difference approximation
!
! pars(npars)          :: input, real values, the parameters of the function
! ED (nED)             :: input, real values, the log-scale OSL ED values 
! error(nED)           :: input, real values, the OSL ED's absolute error
! npars                :: input, integer, the size of ED data
! nED                  :: input, integer, the dimension of the parameters 
! tol                  :: input, real value, tolerance value for diagnosing sigular matrix
! minAbspar            :: input, real value, the allowed minimum absolute parameter 
! hessian(npars,npars) :: output, real values, the hessian matrix
! gradient(npars)      :: output, real values, the gradient of the parameters
! value                :: output, real value, the correspond function value for the specified parameters 
! errorflag(3)         :: output, integer values, if value can be calculated, errorflag(1)=0, otherwise 1;
!                         if gradient can be calculated, errorflag(2)=0, otherwise 1; if hessian can be 
!                         calculated, errorflag(3)=0, otherwise 1. but if sigular matrix appears when 
!                         attempt to call subroutine GJordan to approximate value, gradient and hessian, 
!                         all values in errorflag will be 1
!
! Dependency:: subroutine GJordan, inter function fun34
!
! Author:: Peng Jun, 2013.01.27, revised in 2013.03.17
!
! Reference:  Jose Pinheiro, Douglas Bates, Saikat DebRoy, Deepayan Sarkar and the
!             R Development Core Team (2013). nlme: Linear and Nonlinear Mixed
!             Effects Models. R package version 3.1-108.

!
!--------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::npars                 
  integer(kind=4),intent(in)::nED                    
  integer(kind=4),intent(out)::errorflag(3)             
  real   (kind=8),intent(in)::pars(npars)            
  real   (kind=8),intent(in)::ED(nED)                
  real   (kind=8),intent(in)::Error(nED)            
  real   (kind=8),intent(in)::tol                   
  real   (kind=8),intent(in)::minAbsPar              
  real   (kind=8),intent(out)::gradient(npars)       
  real   (kind=8),intent(out)::hessian(npars,npars)  
  real   (kind=8),intent(out)::value                
  ! local variables
  integer(kind=4)::i,j,p
  integer(kind=4)::ncols
  integer(kind=4)::solerror
  real  (kind=8)::incr(npars)
  real  (kind=8)::diagpar(npars,npars)
  real  (kind=8),allocatable::frac(:),cfrac(:)
  real  (kind=8)::ffrac(1+2*npars)
  real  (kind=8),allocatable::cols(:,:),ccols(:,:),shifted(:,:)
  real  (kind=8)::pcols(npars,2*npars+1)
  real  (kind=8),allocatable::xcols(:,:),cxcols(:,:)
  real  (kind=8),allocatable::pxcols(:,:)
  real(kind=8),parameter::eps=2.013409D-05 ![.Machine$double.eps^0.3 in R]                    
  !
  errorflag=0
  !
  do i=1,npars
    if( abs(pars(i))<=minAbsPar )  then
      incr(i)=minAbsPar*eps
    else
      incr(i)=abs(pars(i))*eps
    end if
  end do
  !
  diagpar=0.0D+00
  do i=1,npars
    diagpar(i,i)=1.0D+00
  end do
  !
  ffrac(1)=1.0D+00
  ffrac(2:npars+1)=incr
  ffrac(npars+2:2*npars+1)=incr**2
  !
  pcols(:,1)=0.0D+00
  pcols(:,2:(npars+1))=diagpar
  pcols(:,(npars+2):(2*npars+1))=-diagpar
  !
  allocate(cols(1:npars,1:(npars*(npars-1)/2)))
  !
  allocate(frac(1:(npars*(npars-1)/2)))
  !
  ncols=0
  do i=1,npars-1
    !
    allocate(ccols(1:npars,1:npars-i))
    allocate(cfrac(1:npars-i))
    !
    do j=i+1,npars
      ccols(:,j-i)=diagpar(:,i)+diagpar(:,j)
      cfrac(j-i)=incr(i)*incr(j)       
    end do 
    !
    cols(:,ncols+1:ncols+npars-i)=ccols
    frac(ncols+1:ncols+npars-i)=cfrac
    !
    deallocate(ccols)
    deallocate(cfrac)
    !
    ncols=ncols+npars-i 
    !
  end do
  !
  allocate(ccols(1:npars,1:npars*(npars-1)/2+2*npars+1))
  ccols(:,1:(2*npars+1))=pcols
  ccols(:,2*npars+2:npars*(npars-1)/2+2*npars+1)=cols
  deallocate(cols)
  allocate(cols(1:npars,1:npars*(npars-1)/2+2*npars+1))
  cols=ccols
  deallocate(ccols)
  !
  allocate(cfrac(1:npars*(npars-1)/2+2*npars+1))
  cfrac(1:2*npars+1)=ffrac
  cfrac(2*npars+2:npars*(npars-1)/2+2*npars+1)=frac
  deallocate(frac)
  allocate(frac(1:npars*(npars-1)/2+2*npars+1))
  frac=cfrac
  deallocate(cfrac) 
  ! 
  p=size(cols,dim=2)
  allocate(shifted(1:npars,p))
  do i=1,p
    shifted(:,i)=pars+incr*cols(:,i)
  end do
  !
  cols=transpose(cols)
  !
  allocate(pxcols(p,1+2*npars))
  pxcols(:,1)=1
  pxcols(:,2:npars+1)=cols
  pxcols(:,npars+2:2*npars+1)=cols**2
  !
  ncols=0
  !
  allocate(xcols(p,1:npars*(npars-1)/2))
  !
  do i=1,npars-1
    allocate(cxcols(p,1:npars-i))
    do j=i+1,npars
      cxcols(:,j-i)=cols(:,i)*cols(:,j)
    end do
    xcols(:,ncols+1:ncols+npars-i)=cxcols
    deallocate(cxcols)
    !
    ncols=ncols+npars-i  
  end do
  !
  allocate(cxcols(p,1:1+2*npars+npars*(npars-1)/2))
  cxcols(:,1:1+2*npars)=pxcols
  cxcols(:,2+2*npars:1+2*npars+npars*(npars-1)/2)=xcols
  deallocate(xcols)
  deallocate(pxcols)
  allocate(xcols(p,1:1+2*npars+npars*(npars-1)/2))
  xcols=cxcols
  deallocate(cxcols)
  !
  allocate(pxcols(1:p,1))
  !
  do i=1,p
    pxcols(i,:)=fun34(shifted(:,i))
  end do
  !
  deallocate(shifted)
  !
  call GJordan(xcols,pxcols,p,1,solerror,tol)
  if(solerror==1)  errorflag=1
  !
  pxcols(:,1)=pxcols(:,1)/frac
  deallocate(xcols)
  deallocate(frac)
  ncols=2*npars+2 
  !
  do j=1,npars
    do i=1,npars
        if (i==j) diagpar(i,j)=pxcols(1+npars+i,1)
        if (i>j ) diagpar(j+1:npars,j)=pxcols(ncols:ncols+npars-j-1,1)  
    end do
    ncols=ncols+npars-j
  end do
  !
  do i=1,npars
    gradient(i)=pxcols(1+i,1)
  end do
  !
  value=pxcols(1,1)
  !
  deallocate(pxcols)
  !
  hessian=diagpar+transpose(diagpar)
  !
  ! check Inf and NaN for value
  if( value .ne. value .or. &
      value+1.0D+00==value  )            errorflag(1)=1
  ! check Inf and NaN for gradient
  if( any(gradient .ne. gradient) .or. &
      any(gradient+1.0D+00==gradient) )  errorflag(2)=1
  ! check Inf and NaN for hessian
  if( any(hessian .ne. hessian) .or. &
      any(hessian+1.0D+00==hessian) )    errorflag(3)=1
  !
  return
  !
  contains
!----------------------
  function fun34(x)
  !------------------------------------------------------------------------------------------------------------
  ! fun34 is a inner function contained in subroutine hessian its used for calculating minus 
  ! logged maximum likelihood value of minimum age model of three or four parameters 
  !
  ! x(npars), input     :: real values, the parameters used for calculating
  ! fun34,   output     :: real value, the value for the specified function
  !
  !  Author :: Peng Jun, 2013.03.16
  !
  ! Reference:: Galbraith, R.F., Roberts, R.G., Laslett, G.M., Yoshida, H. & Olley, J.M., 1999. Optical dating of
  !             single grains of quartz from Jinmium rock shelter, northern Australia. Part I: experimental design
  !             and statistical models. Archaeometry, 41, pp. 339-364.
  !
  ! Dependence:: subroutine pnorm, function alnorm, also the ED data (log-scale) provided in subroutine hessian
  !
  !--------------------------------------------------------------------------------------------------------------
    implicit none
    real(kind=8)::x(npars)
    real(kind=8)::fun34
    real(kind=8)::pnorm1(nED)
    real(kind=8)::alnorm
    real(kind=8),parameter::pi=3.141592653589793238462643383279502884197D+00
    logical,parameter::upper=.false.   
    !
    if (npars==3)    then
      !
      pnorm1=(x(2)-(x(2)/x(3)**2+ED/Error**2)/(1.0D+00/x(3)**2+1.0D+00/Error**2))*sqrt(1.0D+00/Error**2+1.0D+00/x(3)**2)
      !
      call pnorm(pnorm1,nED,upper)
      !
      fun34=-sum(log(x(1)/sqrt(2.0D+00*pi*Error**2)*exp(-(ED-x(2))**2/(2.0D+00*Error**2))+&
	        (1.0D+00-x(1))/sqrt(2.0D+00*pi*(Error**2+x(3)**2))*exp(-(ED-x(2))**2/(2.0D+00*(Error**2+x(3)**2)))*&
	        (1.0D+00-pnorm1)/(0.5D+00)))  
    elseif(npars==4)  then
      !
      pnorm1=(x(2)-(x(3)/x(4)**2+ED/Error**2)/(1.0D+00/x(4)**2+1.0D+00/Error**2))*sqrt(1.0D+00/Error**2+1.0D+00/x(4)**2)
      !
      call pnorm(pnorm1,nED,upper)	 
      !
      fun34=-sum(log(x(1)/sqrt(2.0D+00*pi*Error**2)*exp(-(ED-x(2))**2/(2.0D+00*Error**2))+&
	        (1.0D+00-x(1))/sqrt(2.0D+00*pi*(Error**2+x(4)**2))*exp(-(ED-x(3))**2/(2.0D+00*(Error**2+x(4)**2)))*&
	        (1.0D+00-pnorm1)/(1.0D+00-alnorm((x(2)-x(3))/x(4),upper))))   
      ! 
    end if
    !
    return
    !
  end function fun34   
!---------------------------
end subroutine fdhessian
