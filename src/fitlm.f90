subroutine fitlm(ncomp,tim,sig,ntim,pars,Stdpars,&
                 value,predtval,errorflag)
!---------------------------------------------------------------------------------------------------------------------
! subroutine fitlm is used for fitting OSL signal of type 'lm'
! using Levenberg-Marquadt method, a series combination of 
! initial parameters will be given to call subroutine lmfit
!
! ncomp,             input:: integer, number of components
! ntim,              input:: integer, length of tim or sig
! tim(ntim),         input:: real values, time series
! sig(ntim),         input:: real values, signal series
! pars(2*ncomp),    output:: real values, estimated parameters
! Stdpars(2*ncomp), output:: real values, estimated standard errors for parameters
! value,            output:: real value, minimized objective function value
! predtval(ntim),   output:: real values, predicted signal values correspond to sig
! errorfalg,        output:: integer, error message generated during the calculation, 0 for successful work, 
!                            1 for totally fail, -1 for Std.pars are not available
!
!     Author:: Peng Jun, 2013.07.24, revised in 2013.08.03
!
! Dependence:: subroutine comb_next; subroutine targfunc; subroutine lmfit
! 
! References:: Bluszcz, A., 1996. Exponential function fitting to TL growth data 
!              and similar applications. Geochronometria 13, 135–141.
!              
!              Bulur, E., 2000. A simple transformation for converting CW-OSL 
!              curves to LM-OSL curves. Radiation Measurements 32, 141-145.
!
!              Jain, M., Murray, A.S., Boetter-Jensen, L., 2003. Characterisation of blue-light stimulated 
!              luminescence components in different quartz samples: implications for dose measurement. Radiation
!              Measurements, 37 (4-5), pp. 441-449.
!---------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(kind=4),intent(in)::ncomp
  integer(kind=4),intent(in)::ntim 
  real   (kind=8),dimension(ntim),intent(in)::tim
  real   (kind=8),dimension(ntim),intent(in)::sig
  real   (kind=8),dimension(ntim),intent(out)::predtval
  real   (kind=8),dimension(2*ncomp),intent(out)::pars
  real   (kind=8),dimension(2*ncomp),intent(out)::Stdpars
  integer(kind=4),intent(out)::errorflag
  real   (kind=8),intent(out)::value
  ! variables for subroutine lmfit
  real   (kind=8),dimension(2*ncomp)::lmpars,cpars
  real   (kind=8),dimension(2*ncomp)::lmStdpars,cStdpars
  real   (kind=8),dimension(ntim)::lmpredtval,cpredtval
  real   (kind=8),parameter::lmtol=1.0D-07
  real   (kind=8)::lmvalue,cvalue
  integer(kind=4)::lmErr
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
  integer(kind=4), parameter::typ=2
  integer(kind=4):: i
  !
  errorflag=1
  done=.true.
  Loop: do i=1,permdex(ncomp)
    ! obtain random set
    call comb_next(done,7,ncomp,iarray)
    !
    tlamda=initry(iarray)
    ! obtain tithn
    call targfunc(tlamda,ncomp,tim,sig,targtol,&
                  typ,ntim,tvalue,tithn,targErr)
    ! call lmfit if targErr=0
    if(targErr==0) then
      lmpars=(/tithn,tlamda/)
      call lmfit(tim,sig,ntim,lmpars,2*ncomp,typ,&
                 lmStdpars,lmpredtval,lmvalue,lmtol,lmErr)
      ! reset pars, stdpars, predtval and value if lmErr=1
      if(lmErr==1) then
        pars=lmpars
        Stdpars=lmStdpars
        predtval=lmpredtval
        value=lmvalue
        ! successful simple trails 
        errorflag=0
      end if
      ! if parameters' standard errors can not be estimated, save it too
      if(errorflag/=0 .and. (lmErr==4 .or. lmErr==5 .or. lmErr==6)) then
        cpars=lmpars
        cStdpars=lmStdpars
        cpredtval=lmpredtval
        cvalue=lmvalue
        errorflag=-1
      end if
      ! exit if a successful simple tries has been achieved
      if(errorflag==0) exit Loop
    end if
  end do Loop
  ! if errorflag=0, then return values will be areadly 
  ! saved, else return values need to be reset
  if(errorflag/=0) then
    if(errorflag==-1) then
      ! return a situation that Std.pars can not obtained
      pars=cpars
      Stdpars=cStdpars
      predtval=cpredtval
      value=cvalue
    else if(errorflag==1) then
      ! return a situation that is totally not possible
      pars=    -99.0D+00
      Stdpars= -99.0D+00
      predtval=-99.0D+00
      value=   -99.0D+00
    end if
  end if
  ! now return
  return
end subroutine fitlm
