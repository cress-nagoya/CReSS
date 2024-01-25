!***********************************************************************
      module m_setvar1d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/10,
!                   1999/05/20, 1999/06/21, 1999/06/28, 1999/07/05,
!                   1999/07/23, 1999/08/03, 1999/08/18, 1999/08/23,
!                   1999/09/06, 1999/09/30, 1999/10/12, 1999/11/01,
!                   1999/12/06, 2000/01/17, 2000/03/17, 2000/12/18,
!                   2001/01/15, 2001/03/13, 2001/04/15, 2001/06/06,
!                   2001/07/13, 2001/08/07, 2001/09/13, 2002/04/02,
!                   2002/06/06, 2002/08/15, 2002/10/31, 2003/03/21,
!                   2003/05/19, 2003/12/12, 2004/05/31, 2004/08/20,
!                   2004/10/12, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/07/21, 2007/01/20, 2007/01/31, 2007/05/14,
!                   2007/07/30, 2007/11/26, 2008/05/02, 2008/08/25,
!                   2009/01/30, 2009/02/27, 2011/05/16, 2011/08/18,
!                   2011/09/22, 2011/11/10, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the initial velocity and the scalar perturbations.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_commath
      use m_copy3d
      use m_getcname
      use m_getiname
      use m_getpt0
      use m_getqt0
      use m_inichar
      use m_setcst3d
      use m_setcst4d
      use m_shift2nd

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setvar1d, s_setvar1d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setvar1d

        module procedure s_setvar1d

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_setvar1d(fpdmpvar,fpadvopt,fpcphopt,fpqcgopt,        &
     &                      fpaslopt,fptrkopt,fptubopt,fmois,           &
     &                      ni,nj,nk,nqw,nnw,nqi,nni,nqa,zph,ubr,vbr,   &
     &                      qvbr,u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,qv,qvp, &
     &                      qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,&
     &                      qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,&
     &                      tke,tkep,maxvl,prwtr,price)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpqcgopt
                       ! Formal parameter of unique index of qcgopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

      integer, intent(in) :: fptrkopt
                       ! Formal parameter of unique index of trkopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nni
                       ! Number of categories of ice concentrations

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

! Input and output variables

      real, intent(inout) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(inout) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(inout) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

! Output variables

      real, intent(out) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(out) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(out) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(out) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(out) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(out) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(out) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at present

      real, intent(out) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(out) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(out) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(out) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(out) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(out) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at present

      real, intent(out) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(out) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at present

      real, intent(out) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(out) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

      real, intent(out) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(out) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at present

      real, intent(out) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(out) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at present

      real, intent(out) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(out) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at present

      real, intent(out) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(out) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(out) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(out) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at present

      real, intent(out) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

      real, intent(out) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(out) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(out) :: maxvl(0:ni+1,0:nj+1,1:nk)
                       ! Maximum instantaneous wind velocity

      real, intent(out) :: prwtr(0:ni+1,0:nj+1,1:2,1:nqw)
                       ! Precipitation and accumulation for water

      real, intent(out) :: price(0:ni+1,0:nj+1,1:2,1:nqi)
                       ! Precipitation and accumulation for ice

! Internal shared variables

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      integer n        ! Array index in 4th direction

! Remark

!     up,vp: These variables are also temporary.

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(dmpvar)

! -----

! Get the required namelist variables.

      call getcname(fpdmpvar,dmpvar)
      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)

! -----

! Set the initial potential temperature perturbation.

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ptpp)

      call s_getpt0(idpt0opt,idpt0num,idptp0,idpt0rx,idpt0ry,idpt0rz,   &
     &              idpt0cx,idpt0cy,idpt0cz,idpt0ds,ni,nj,nk,           &
     &              zph,ptpp,up,vp)

! -----

! Set the initial aerosol.

      if(aslopt.ge.1) then

        do n=1,nqa(0)

!ORIG     call s_setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qaslp(0,0,1,n))

          call s_getqt0(idqt0opt,idqt0num,idqt0,idqt0rx,idqt0ry,idqt0rz,&
     &                  idqt0cx,idqt0cy,idqt0cz,idqt0ds,ni,nj,nk,zph,   &
     &                  qaslp(0,0,1,n),up,vp)

        end do

      end if

! -----

! Set the initial tracer.

      if(trkopt.ge.1) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qtp)

        if(trkopt.eq.1) then

          call s_getqt0(idqt0opt,idqt0num,idqt0,idqt0rx,idqt0ry,idqt0rz,&
     &                  idqt0cx,idqt0cy,idqt0cz,idqt0ds,ni,nj,nk,zph,   &
     &                  qtp,up,vp)

        end if

      end if

! -----

! Exchange the value in the case the 4th order calculation is performed.

      call shift2nd(idwbc,idebc,idsbc,idnbc,idadvopt,idsmtopt,idcphopt, &
     &              idhaiopt,idqcgopt,idaslopt,idtrkopt,idtubopt,fmois, &
     &              'ooxxooxoox',ni,nj,nk,nqw,nnw,nqi,nni,nqa,ubr,vbr,  &
     &              wp,ppp,ptpp,qvbr,qwtrp,nwtrp,qicep,nicep,           &
     &              qcwtrp,qcicep,qaslp,qtp,tkep)

! -----

! Finally set all of the prognostic variables.

      if(advopt.le.3) then

        call copy3d(0,ni+1,0,nj+1,1,nk,ubr,u)
        call copy3d(0,ni+1,0,nj+1,1,nk,ubr,up)

        call copy3d(0,ni+1,0,nj+1,1,nk,vbr,v)
        call copy3d(0,ni+1,0,nj+1,1,nk,vbr,vp)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,w)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,wp)

        call copy3d(0,ni+1,0,nj+1,1,nk,ptpp,ptp)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,pp)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ppp)

        call copy3d(0,ni+1,0,nj+1,1,nk,qvbr,qv)
        call copy3d(0,ni+1,0,nj+1,1,nk,qvbr,qvp)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwtr)
            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwtrp)

            if(abs(cphopt).eq.4) then

              call setcst4d(0,ni+1,0,nj+1,1,nk,1,nnw,0.e0,nwtr)

            end if

            if(abs(cphopt).ge.2) then

              call setcst4d(0,ni+1,0,nj+1,1,nk,1,nnw,0.e0,nwtrp)

            end if

            call setcst4d(0,ni+1,0,nj+1,1,2,1,nqw,0.e0,prwtr)

          end if

          if(abs(cphopt).ge.2) then

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qice)
            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qicep)

            if(abs(cphopt).eq.2) then

              call setcst4d(0,ni+1,0,nj+1,1,nk,1,1,0.e0,nice)

            else

              call setcst4d(0,ni+1,0,nj+1,1,nk,1,nni,0.e0,nice)

            end if

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nni,0.e0,nicep)

            call setcst4d(0,ni+1,0,nj+1,1,2,1,nqi,0.e0,price)

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qcwtr)
              call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qcwtrp)

            end if

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qcice)
            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qcicep)

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwtr)
            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwtrp)

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nnw,0.e0,nwtr)
            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nnw,0.e0,nwtrp)

            call setcst4d(0,ni+1,0,nj+1,1,2,1,1,0.e0,prwtr)

          end if

          if(abs(cphopt).eq.12) then

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qice)
            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qicep)

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nni,0.e0,nice)
            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nni,0.e0,nicep)

            call setcst4d(0,ni+1,0,nj+1,1,2,1,1,0.e0,price)

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            call s_copy3d(0,ni+1,0,nj+1,1,nk,qaslp(0,0,1,n),            &
     &                    qasl(0,0,1,n))

          end do

        end if

        if(trkopt.ge.1) then

          call copy3d(0,ni+1,0,nj+1,1,nk,qtp,qt)

        end if

        if(tubopt.ge.2) then

          call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tke)
          call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tkep)

          if(dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-') then

            call setcst3d(0,ni+1,0,nj+1,1,nk,lim35n,maxvl)

          end if

        end if

      else

        call copy3d(0,ni+1,0,nj+1,1,nk,ubr,up)
        call copy3d(0,ni+1,0,nj+1,1,nk,vbr,vp)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,wp)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ppp)

        call copy3d(0,ni+1,0,nj+1,1,nk,qvbr,qvp)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwtrp)

            if(abs(cphopt).ge.2) then

              call setcst4d(0,ni+1,0,nj+1,1,nk,1,nnw,0.e0,nwtrp)

            end if

            call setcst4d(0,ni+1,0,nj+1,1,2,1,nqw,0.e0,prwtr)

          end if

          if(abs(cphopt).ge.2) then

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qicep)
            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nni,0.e0,nicep)

            call setcst4d(0,ni+1,0,nj+1,1,2,1,nqi,0.e0,price)

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qcwtrp)

            end if

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qcicep)

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwtrp)
            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nnw,0.e0,nwtrp)

            call setcst4d(0,ni+1,0,nj+1,1,2,1,1,0.e0,prwtr)

          end if

          if(abs(cphopt).eq.12) then

            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qicep)
            call setcst4d(0,ni+1,0,nj+1,1,nk,1,nni,0.e0,nicep)

            call setcst4d(0,ni+1,0,nj+1,1,2,1,1,0.e0,price)

          end if

        end if

        if(tubopt.ge.2) then

          call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tkep)

          if(dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-') then

            call setcst3d(0,ni+1,0,nj+1,1,nk,lim35n,maxvl)

          end if

        end if

      end if

! -----

      end subroutine s_setvar1d

!-----7--------------------------------------------------------------7--

      end module m_setvar1d
