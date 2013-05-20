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
! Author:: Peng Jun, 2013.01.27, revised in 2013.03.17, revised again in 2013.04.21
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
  real  (kind=8),allocatable::cols(:,:),ccols(:,:),shifted(:,:),transcols(:,:)
  real  (kind=8)::pcols(npars,2*npars+1)
  real  (kind=8),allocatable::xcols(:,:),cxcols(:,:)
  real  (kind=8),allocatable::pxcols(:,:)
  real(kind=8),parameter::eps=2.013409D-05 ![.Machine$double.eps^0.3 in R]                    
  !
  errorflag=0
  ! Decide the incr values
  ! for each initial par in pars, check if it
  ! is smaller than minabspar, the incr will
  ! be decided through this check
  do i=1,npars
    if( abs(pars(i))<=minAbsPar )  then
      incr(i)=minAbsPar*eps
    else
      incr(i)=abs(pars(i))*eps
    end if
  end do
  !
  ! build a diagal matirx diagpar
  diagpar=0.0D+00
  do i=1,npars
    diagpar(i,i)=1.0D+00
  end do
  ! creat ffrac with length of 2*npars+1
  ! and specify values for it
  ffrac(1)=1.0D+00
  ffrac(2:npars+1)=incr
  ffrac(npars+2:2*npars+1)=incr**2
  ! create pcols to be npars rows
  ! and 2*npars columns, then storing
  ! diagpar in it
  pcols(:,1)=0.0D+00
  pcols(:,2:npars+1)=diagpar
  pcols(:,npars+2:2*npars+1)=-diagpar
  ! in a total of (npars-1) times looping
  ! add  (npars-i) columns to cols in each loop number
  ! hence a toal of npars*(npars-1)/2 columns for cols
  allocate( cols(1:npars,1:npars*(npars-1)/2) )
  ! in a total of (npars-1) times looping
  ! add (npars-i) number to frac in each loop number
  ! hence a total of npars*(npars-1)/2 values for frac
  allocate( frac(1:npars*(npars-1)/2) )
  !
  ncols=0
  do i=1,npars-1
    ! in each loop(i), generate ccols and cfrac with 
    ! npars rows (npars-i) columns 
    allocate( ccols(1:npars,1:npars-i) )
    allocate( cfrac(1:npars-i) )
    ! fill ccols and cfrac in each loop
    do j=i+1,npars
      ccols(:,j-i)=diagpar(:,i)+diagpar(:,j)
      cfrac(j-i)=incr(i)*incr(j)       
    end do 
    ! now store ccols and cfrac to cols
    ! and frac respectively
    cols(:,ncols+1:ncols+npars-i)=ccols
    frac(ncols+1:ncols+npars-i)=cfrac
    ! delocate ccols and cfrac, preparing
    ! for new allocations
    deallocate(ccols)
    deallocate(cfrac)
    ! update started filling index ncols
    ncols=ncols+npars-i 
    !
  end do
  ! now add up pcols and cols togeteher to be ccols
  ! ccols has a total of ( 2*npars+1 + npars*(npars-1)/2 ) columns
  allocate( ccols(1:npars,1:npars*(npars-1)/2+2*npars+1) )
  ccols(:,1:2*npars+1)=pcols
  ccols(:,2*npars+2:npars*(npars-1)/2+2*npars+1)=cols
  ! not need cols presently, so release it
  deallocate(cols)
  ! allocate cols again to store ccols, and release ccols
  allocate(cols(1:npars,1:npars*(npars-1)/2+2*npars+1))
  cols=ccols
  deallocate(ccols)
  ! now add up ffrac and frac together to be cfrac,
  ! cfrac has a total of ( 2*npars+1 + npars*(npars-1)/2 ) values
  allocate( cfrac(1:npars*(npars-1)/2+2*npars+1) )
  cfrac(1:2*npars+1)=ffrac
  cfrac(2*npars+2:npars*(npars-1)/2+2*npars+1)=frac
  ! not need frac so release it
  deallocate(frac)
  ! allocate frac again, store cfrac in it then release cfrac
  allocate(frac(1:npars*(npars-1)/2+2*npars+1))
  frac=cfrac
  deallocate(cfrac) 
  ! specify p to be the column number of cols
  ! that is npars*(npars-1)/2+2*npars+1
  p=npars*(npars-1)/2+2*npars+1
  ! now allocate shifted to be the same shape 
  ! with cols and store some values in it
  allocate(shifted(1:npars,p))
  do i=1,p
    shifted(:,i)=pars+incr*cols(:,i)
  end do
  ! allocate transcols to store transposed cols
  ! transcols has p rows, npars columns
  allocate( transcols(1:p,1:npars) )
  transcols=transpose(cols)
  ! now cols will not be needed, release it
  deallocate(cols)
  ! allocate pxcols to store transcols 
  ! in differ style
  allocate(pxcols(p,1+2*npars))
  pxcols(:,1)=1.0D+00
  pxcols(:,2:npars+1)=transcols
  pxcols(:,npars+2:2*npars+1)=(transcols)**2
  !
  ncols=0
  ! now allocate xcols with p rows and 
  ! npars*(npars-1)/2 columns
  allocate( xcols(p,1:npars*(npars-1)/2) )
  ! 
  do i=1,npars-1
    ! in each loop, allocate cxcols 
    allocate(cxcols(p,1:npars-i))
    do j=i+1,npars
      cxcols(:,j-i)=transcols(:,i)*transcols(:,j)
    end do
    ! store cxcols to xcols and release it
    xcols(:,ncols+1:ncols+npars-i)=cxcols
    deallocate(cxcols)
    ! update started filling index
    ncols=ncols+npars-i  
  end do
  ! now store pxcols and xcols together to cxcols
  allocate( cxcols(p,1:1+2*npars+npars*(npars-1)/2) )
  cxcols(:,1:1+2*npars)=pxcols
  cxcols(:,2+2*npars:1+2*npars+npars*(npars-1)/2)=xcols
  ! release xcols and pxcols
  deallocate(xcols)
  deallocate(pxcols)
  ! allocate xcols again to store cxcols
  allocate(xcols(p,1:1+2*npars+npars*(npars-1)/2))
  xcols=cxcols
  ! release cxcols
  deallocate(cxcols)
  ! allocate pxcols again and store some
  ! values to it then release shifted
  ! now pxcols has p rows and 1 column
  allocate(pxcols(1:p,1))
  do i=1,p
    pxcols(i,:)=fun34(shifted(:,i))
  end do
  deallocate(shifted)
  ! solve xcols %*% X = pxcols
  ! store solved X values in pxcols
  call GJordan(xcols,pxcols,p,1,solerror,tol)
  if(solerror==1)  errorflag=1
  ! scale pxcols with frac and release frac,
  ! release xcols
  pxcols(:,1)=pxcols(:,1)/frac
  deallocate(xcols)
  deallocate(frac)
  ncols=2*npars+2 
  ! fill diagpar with some new values,  note that
  ! non-diagnal part of digpar are zeros
  ! these value are comes from pxcols, and will
  ! be used to constructe hessian matrix
  do j=1,npars
    do i=1,npars
        if (i==j) diagpar(i,j)=pxcols(1+npars+i,1)
        if (i>j ) diagpar(j+1:npars,j)=pxcols(ncols:ncols+npars-j-1,1)  
    end do
    ncols=ncols+npars-j
  end do
  ! estimate gradients
  do i=1,npars
    gradient(i)=pxcols(1+i,1)
  end do
  ! estimate fun(pars)
  value=pxcols(1,1)
  ! now pxcols will not be needed, release it
  deallocate(pxcols)
  ! estimate hessian matrix
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
  ! if no error appears, return 
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
