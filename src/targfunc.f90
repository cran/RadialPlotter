subroutine targfunc(lamda,nlamda,tim,sig,tol,&
                    ntim,value,ithn,errorflag)
!---------------------------------------------------------------------------------------------------------
! targfunc is an objective function for hela using;
! initial lamda values must be provided, then ithn values
! will be calculated using Linear Algebra method
! I(t)=ithn1*exp(-lamda1*t)+ithn2*exp(-lamda2*t)+...+ithnk*exp(-lamdak*t)
!
! lamda(nlamda),     input:: real values, lamda values
! nlamda,            input:: integer, length of lamda values
! tim(ntim),         input:: real values, time values
! sig(ntim),         input:: real values, signal values
! tol,               input:: real value, maximum tolerance for lenear independence
! ntim,              input:: integer, length of tim and sig
! value,            output:: function value
! ithn(nlamda),     output:: I0i values (i=1,2,...maxcomp)
! errorflag,        output:: integer, error message generated during the calculation
!
! Author:: Peng Jun, 2013.05.20, revised in 2013.05.31
!
! Dependence:: subroutine GJordan
!
! References:: Bluszcz, A., 1996. Exponential function fitting to TL growth data and similar
!              applications. Geochronometria 13, 135–141.
!
!              Bluszcz, A., Adamiec, G., 2006. Application of differential evolution to fitting
!              OSL decay curves. Radiation Measurements 41, 886-891.
!----------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::nlamda                        ! length of lamda values
  integer(kind=4),intent(in)::ntim                          ! length of tim and sig
  integer(kind=4),intent(out)::errorflag                    ! error message generaged during the calculation
  real   (kind=8),intent(in)::tol                           ! maximum tolerance for linear independence
  real   (kind=8),dimension(nlamda),intent(in)::lamda       ! lamda values
  real   (kind=8),dimension(ntim),intent(in)::tim,sig       ! time and signal values
  real   (kind=8),intent(out)::value                        ! targeted function value
  real   (kind=8),dimension(nlamda),intent(out)::ithn       ! I0i values (i=1,2,...maxcomp)
  ! local variables
  real(kind=8),dimension(ntim,nlamda)::coef                 ! coefficent matrix 
  real(kind=8),dimension(ntim,1)::ssig                      ! for storing sig values using matrix
  real(kind=8),dimension(nlamda,1)::iithn                   ! for storing ithn values using matrix
  real(kind=8),dimension(ntim)::rowsumcoef                  ! sum values in matrix coef by row
  real(kind=8),dimension(nlamda,nlamda)::ccoef              ! t(coef)%*%coef 
  integer(kind=4)::i                                        ! iterative lopping index
  !
  ! initializing coefficents 
  do i=1,nlamda
    coef(:,i)=dexp( -lamda(i)*tim(:) )
  end do
  !
  ! store sig values in matrix ssig
  ssig(:,1)=sig
  !
  ! calculate ithn values using Linear Algebra method
  ! initialize iithn and ccoef for Linear Algebra using
  iithn=matmul(transpose(coef),ssig)
  ccoef=matmul(transpose(coef),coef)
  call GJordan(ccoef,iithn,nlamda,1,errorflag,tol)
  ! at this point, if no error appears when calling
  ! GJordan, errorflag will be 0
  !
  ! pass values of iithn to ithn
  ithn=iithn(:,1)
  !
  ! return if error appears. error situations including:
  ! 1) matrix ccoef is singular; 2) any obtained ithn value is below zero
  ! then errorflag=1, ithn=ithn(wrong values), value=1.00D+30
  if(errorflag==1 .or. any(ithn<1.0D-10) )  then
    value=1.00D+30
    errorflag=1
    return
  end if
  !
  ! reset matirx coef
  do i=1,nlamda
    coef(:,i)=ithn(i)*dexp( -lamda(i)*tim )
  end do
  !
  ! calculate rowsums of matrix coef
  do i=1,ntim
    rowsumcoef(i)=sum(coef(i,:))
  end do
  ! set residual value
  value=sum( (sig-rowsumcoef)**2 ) 
  return
end subroutine targfunc
