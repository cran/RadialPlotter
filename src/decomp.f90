subroutine decomp(ncomp,tim,sig,ntim,pars,Stdpars,value,&
                  predtval,factor,f,cr,maxiter,tol,errorflag)
!----------------------------------------------------------------------------------------------------------------------------
! subroutine decomp is used to decompose CW OSL decay curve
!
! ncomp,             input:: integer, the number of components to be decomposed
! tim(ntim),         input:: real values, time values
! sig(ntim),         input:: real values, signal values
! ntim,              input:: integer, length of signal values
! pars(2*ncomp),    output:: real values, calculated parameters
! Stdpars(2*ncomp), output:: real values, calculated standard errors for parameters
! value,            output:: real value, minimized sum of square of residuals
! predtval(ntim),   output:: real values, fitted values correspond to sig values
! factor,            input:: integer, for scaling np in differential evolution, np=ncomp*factor
! f   ,              input:: real value, parameter for differential evolution
! cr,                input:: real value, parameter for differential evolution
! maxiter,           input:: integer, maximum allowed iteration for differential evolution
! tol,               input:: real value, tolerance for stoping iteration in differential evolution
! errorflag(3),     output:: integer values, error message generated during the running:
!                            1) if parameters can not be initialized using differential evolution, errorflag(1)=1, else 0;
!                            2) if parameters can not be optimized using Levenberg-Marquadt method, errorflag(2)=1, else 0;
!                            3) if simple trails have not been performed, errorflag(3)=0, if success in simple trails, 
!                               errorflag(3)=1,else errorflag(3)=-1.
!
! Author:: Peng Jun, 2013.06.05
!
! Dependence:: subroutine diffev; subroutine lmfit; subroutine comb_next; subroutine targfunc.
!
! References:: Bluszcz, A., 1996. Exponential function fitting to TL growth data and similar
!              applications. Geochronometria 13, 135–141.
!
!              Bluszcz, A., Adamiec, G., 2006. Application of differential evolution to fitting
!              OSL decay curves. Radiation Measurements 41, 886-891.
!     
!              Jain, M., Murray, A.S., Boetter-Jensen, L., 2003. Characterisation of blue-light stimulated 
!              luminescence components in different quartz samples: implications for dose measurement. Radiation
!              Measurements, 37 (4-5), pp. 441-449.
!            
!              Jorge More, Burton Garbow, Kenneth Hillstrom,User Guide for MINPACK-1,
!              Technical Report ANL-80-74, Argonne National Laboratory, 1980.   
!---------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::ncomp
  integer(kind=4),intent(in)::ntim 
  integer(kind=4),intent(in)::factor
  integer(kind=4),intent(in)::maxiter
  real   (kind=8),intent(in)::f
  real   (kind=8),intent(in)::cr
  real   (kind=8),intent(in)::tol
  real   (kind=8),dimension(ntim),intent(in)::tim
  real   (kind=8),dimension(ntim),intent(in)::sig
  real   (kind=8),dimension(ntim),intent(out)::predtval
  real   (kind=8),dimension(2*ncomp),intent(out)::pars
  real   (kind=8),dimension(2*ncomp),intent(out)::Stdpars
  integer(kind=4),dimension(3),intent(out)::errorflag
  real   (kind=8),intent(out)::value
  !
  ! variables for subroutine diffev
  real   (kind=8),dimension(ncomp,factor*ncomp)::agents
  real   (kind=8),dimension(ncomp)::dlamda,dithn
  real   (kind=8),dimension(ntim)::dpredict
  integer(kind=4),dimension(2)::dErr
  real   (kind=8)::dvalue
  !
  ! variables for subroutine lmfit
  real   (kind=8),dimension(2*ncomp)::lmpars
  real   (kind=8),dimension(2*ncomp)::lmStdpars
  real   (kind=8),dimension(ntim)::lmpredict
  real   (kind=8),parameter::lmtol=1.0D-07
  real   (kind=8)::lmvalue
  integer(kind=4)::lmErr
  !
  ! variables for subroutine targfunc and subroutine comb_next
  logical:: done
  real   (kind=8),parameter::targtol=1.0D-07
  integer(kind=4)::targErr
  real   (kind=8)::tvalue
  integer(kind=4),dimension(ncomp)::iarray
  real   (kind=8),dimension(ncomp)::tlamda,tithn
  integer(kind=4),dimension(7),parameter::permdex=(/7,21,35,35,21,7,1/)
  real   (kind=8),dimension(7),parameter::initry=(/32.0D+00,  &      ! UF
                                                   2.5D+00 ,  &      ! Fast 
                                                   0.62D+00,  &      ! Medium 
                                                   0.15D+00,  &      ! Slow1 
                                                   0.023D+00, &      ! Slow2 
                                                   0.0022D+00,&      ! Slow3 
                                                   0.0003D+00 /)     ! Slow4
  ! local variables
  integer(kind=4)::i
  !
  ! at this points all errors are 0
  errorflag=0
  ! 
  ! call diffev
  call diffev(ncomp,tim,sig,ntim,dlamda,dithn,dvalue,dpredict,&
              agents,factor*ncomp,f,cr,maxiter,tol,dErr)
  !
  ! set pars, stdpars, predtval and value if dErr(2)=0
  if(dErr(2)==0) then
    pars=(/dithn,dlamda/)
    Stdpars=-99.0D+00
    predtval=dpredict
    value=dvalue
  else 
    ! set errorflag(1) to be 1 if diffev fails 
    errorflag(1)=1
  end if 
  !
  ! call lmfit if diffev success
  if(errorflag(1)==0) then
    lmpars=(/dithn,dlamda/)
    call lmfit(tim,sig,ntim,lmpars,2*ncomp,&
               lmStdpars,lmpredict,lmvalue,lmtol,lmErr)
    ! reset pars, stdpars, predtval and value if lmErr=0
    if(lmErr==1) then
      pars=lmpars
      Stdpars=lmStdpars
      predtval=lmpredict
      value=lmvalue
    else 
      ! set errorflag(2) to be 1 if lmfit fails
      errorflag(2)=1
    end if
    !
  end if
  !
  ! if diffev or lmfit fails then simple trials will be used
  if(errorflag(1)/=0 .or. errorflag(2)/=0) then
    ! if fail in simple trails, errorflag(3) will be -1
    errorflag(3)=-1
    done=.true.
    Loop: do i=1,permdex(ncomp)
      ! obtain random set
      call comb_next(done,7,ncomp,iarray)
      !
      tlamda=initry(iarray)
      ! obtain tithn
      call targfunc(tlamda,ncomp,tim,sig,targtol,&
                      ntim,tvalue,tithn,targErr)
      ! call lmfit if targErr=0
      if(targErr==0) then
        lmpars=(/tithn,tlamda/)
        call lmfit(tim,sig,ntim,lmpars,2*ncomp,lmStdpars,&
                   lmpredict,lmvalue,lmtol,lmErr)
        ! reset pars, stdpars, predtval and value if lmErr=0
        if(lmErr==1) then
          pars=lmpars
          Stdpars=lmStdpars
          predtval=lmpredict
          value=lmvalue
          ! successful simple trails 
          errorflag(3)=1
        end if
        ! exit if a successful simple trail has been achieved
        if(errorflag(3)==1) exit Loop
      end if
    end do Loop
  end if
  !
  ! if fail both in diffev and simple trails, set output to be -99.0
  if(errorflag(1)==1 .and. errorflag(3)==-1) then
    pars=-99.0D+00
    Stdpars=-99.0D+00
    predtval=-99.0D+00
    value=-99.0D+00
  end if
  !
  return
end subroutine decomp
