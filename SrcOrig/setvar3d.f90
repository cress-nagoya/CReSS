!***********************************************************************
      module m_setvar3d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/06/21, 1999/07/05, 1999/08/03, 1999/08/18,
!                   1999/08/23, 1999/09/06, 1999/10/12, 1999/10/27,
!                   1999/11/01, 1999/11/19, 1999/12/06, 1999/12/20,
!                   2000/01/17, 2000/03/17, 2000/04/18, 2000/06/01,
!                   2000/12/18, 2001/01/15, 2001/03/13, 2001/06/06,
!                   2001/07/13, 2001/08/07, 2001/09/13, 2001/10/18,
!                   2001/11/20, 2002/04/02, 2002/06/06, 2002/08/15,
!                   2002/10/31, 2002/12/02, 2003/01/04, 2003/03/21,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/11/28,
!                   2003/12/12, 2004/03/05, 2004/04/01, 2004/05/31,
!                   2004/06/10, 2004/07/01, 2004/08/20, 2004/09/01,
!                   2004/09/25, 2004/10/12, 2004/12/17, 2005/10/05,
!                   2006/01/10, 2006/02/13, 2006/04/03, 2006/06/21,
!                   2006/07/21, 2006/09/21, 2006/11/06, 2007/01/20,
!                   2007/01/31, 2007/05/07, 2007/05/14, 2007/07/30,
!                   2007/11/26, 2008/05/02, 2008/06/09, 2008/08/25,
!                   2009/01/30, 2009/02/27, 2011/05/16, 2011/08/18,
!                   2011/09/22, 2011/11/10, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the initial and the boundary conditions.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_commath
      use m_copy3d
      use m_copy4d
      use m_diagnci
      use m_diagni
      use m_diagnw
      use m_getcname
      use m_getiname
      use m_getqt0
      use m_inichar
      use m_phy2cnt
      use m_setcst3d
      use m_setcst4d
      use m_shift2nd
      use m_vbcp
      use m_vbcs
      use m_vbcu
      use m_vbcv
      use m_vbcw

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setvar3d, s_setvar3d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setvar3d

        module procedure s_setvar3d

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
      subroutine s_setvar3d(fpgpvvar,fpdmpvar,fpadvopt,fpcphopt,        &
     &                      fphaiopt,fpqcgopt,fpaslopt,fptrkopt,        &
     &                      fptubopt,fmois,ni,nj,nk,nqw,nnw,nqi,nni,    &
     &                      nqa,zph,j31,j32,jcb8w,mf,rbr,up,vp,wp,      &
     &                      ppp,ptpp,qvp,qwtrp,qicep,u,v,w,pp,ptp,      &
     &                      qv,qwtr,nwtr,nwtrp,qice,nice,nicep,         &
     &                      qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,       &
     &                      qt,qtp,tke,tkep,maxvl,prwtr,price,          &
     &                      tmp1,tmp2,tmp3,tmp4)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

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

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

! Input and output variables

      real, intent(inout) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(inout) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(inout) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(inout) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(inout) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(inout) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(inout) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(inout) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

! Output variables

      real, intent(out) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(out) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(out) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(out) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at present

      real, intent(out) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(out) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(out) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at present

      real, intent(out) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at present

      real, intent(out) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(out) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

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

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      integer n        ! Array index in 4th direction

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(gpvvar)
      call inichar(dmpvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getcname(fpdmpvar,dmpvar)
      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)

! -----

! Set the boundary conditions for the x components of velocity.

      call vbcu(ni,nj,nk,up)

! -----

! Set the boundary conditions for the y components of velocity.

      call vbcv(ni,nj,nk,vp)

! -----

! Set the boundary conditions for the z components of velocity.

      if(gpvvar(1:1).eq.'x') then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,wp)

      end if

      call s_phy2cnt(idsthopt,idtrnopt,idmpopt,idmfcopt,idoneopt,       &
     &               ni,nj,nk,j31,j32,jcb8w,mf,up,vp,wp,                &
     &               tmp1,tmp2,tmp3,tmp4)

      call s_vbcw(idbbc,idtbc,idmpopt,idmfcopt,ni,nj,nk,j31,j32,jcb8w,  &
     &            mf,up,vp,tmp1,wp,tmp2,tmp3,tmp4)

! -----

! Set the boundary conditions for the pressure.

      call vbcp(idbbc,ni,nj,nk,ppp)

! -----

! Set the boundary conditions for the potential temperature.

      call vbcs(ni,nj,nk,ptpp)

! -----

! Set the boundary conditions for the water vapor mixing ratio.

      if(gpvvar(2:2).eq.'x') then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qvp)

      else

        call vbcs(ni,nj,nk,qvp)

      end if

! -----

!!! Set the boundary conditions for the water and ice hydrometeor.

!! For the bulk categories.

      if(abs(cphopt).lt.10) then

! Set the boundary conditions for the water hydrometeor.

        if(abs(cphopt).ge.1) then

          if(gpvvar(3:3).eq.'x') then

            call s_setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qwtrp(0,0,1,1))

          else

            call s_vbcs(ni,nj,nk,qwtrp(0,0,1,1))

          end if

          if(gpvvar(4:4).eq.'x') then

            call s_setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qwtrp(0,0,1,2))

          else

            call s_vbcs(ni,nj,nk,qwtrp(0,0,1,2))

          end if

        end if

! -----

! Set the boundary conditions for the ice hydrometeor.

        if(abs(cphopt).ge.2) then

          if(gpvvar(5:5).eq.'x') then

            call s_setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qicep(0,0,1,1))

          else

            call s_vbcs(ni,nj,nk,qicep(0,0,1,1))

          end if

          if(gpvvar(6:6).eq.'x') then

            call s_setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qicep(0,0,1,2))

          else

            call s_vbcs(ni,nj,nk,qicep(0,0,1,2))

          end if

          if(gpvvar(7:7).eq.'x') then

            call s_setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qicep(0,0,1,3))

          else

            call s_vbcs(ni,nj,nk,qicep(0,0,1,3))

          end if

          if(haiopt.eq.1) then

            if(gpvvar(8:8).eq.'x') then

              call s_setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qicep(0,0,1,4))

            else

              call s_vbcs(ni,nj,nk,qicep(0,0,1,4))

            end if

          end if

        end if

! -----

      end if

!! -----

!!! -----

! Calculate the diagnostic concentrations for the cloud ice, snow and
! graupel.

      if(abs(cphopt).lt.10) then

        if(abs(cphopt).eq.2) then

          call diagnci(ni,nj,nk,nqi,nni,qicep,nicep)

        else if(abs(cphopt).eq.3) then

          call diagni(idhaiopt,ni,nj,nk,nqi,nni,rbr,qicep,nicep)

        else if(abs(cphopt).eq.4) then

          call diagnw(ni,nj,nk,nqw,nnw,rbr,qwtrp,nwtrp)

          call diagni(idhaiopt,ni,nj,nk,nqi,nni,rbr,qicep,nicep)

        end if

      end if

! -----

! Set the initial aerosol.

      if(aslopt.ge.1) then

        do n=1,nqa(0)

!ORIG     call s_setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qaslp(0,0,1,n))

!ORIG     call s_getqt0(idqt0opt,idqt0num,idqt0,idqt0rx,idqt0ry,idqt0rz,&
!ORIG&                  idqt0cx,idqt0cy,idqt0cz,idqt0ds,ni,nj,nk,zph,   &
!ORIG&                  qaslp(0,0,1,n),tmp1,tmp2)

          call s_vbcs(ni,nj,nk,qaslp(0,0,1,n))

        end do

      end if

! -----

! Set the initial tracer.

      if(trkopt.ge.1) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qtp)

        if(trkopt.eq.1) then

          call s_getqt0(idqt0opt,idqt0num,idqt0,idqt0rx,idqt0ry,idqt0rz,&
     &                  idqt0cx,idqt0cy,idqt0cz,idqt0ds,ni,nj,nk,zph,   &
     &                  qtp,tmp1,tmp2)

        end if

      end if

! -----

! Exchange the value in the case the 4th order calculation is performed.

      if(abs(cphopt).lt.10) then

       call shift2nd(idwbc,idebc,idsbc,idnbc,idadvopt,idsmtopt,idcphopt,&
     &               idhaiopt,idqcgopt,idaslopt,idtrkopt,idtubopt,fmois,&
     &               'ooooooooox',ni,nj,nk,nqw,nnw,nqi,nni,nqa,up,vp,wp,&
     &               ppp,ptpp,qvp,qwtrp,nwtrp,qicep,nicep,qcwtrp,qcicep,&
     &               qaslp,qtp,tkep)

      else

       call shift2nd(idwbc,idebc,idsbc,idnbc,idadvopt,idsmtopt,idcphopt,&
     &               idhaiopt,idqcgopt,idaslopt,idtrkopt,idtubopt,fmois,&
     &               'ooooooxoox',ni,nj,nk,nqw,nnw,nqi,nni,nqa,up,vp,wp,&
     &               ppp,ptpp,qvp,qwtrp,nwtrp,qicep,nicep,qcwtrp,qcicep,&
     &               qaslp,qtp,tkep)

      end if

! -----

! Finally set all of the prognostic variables.

      if(advopt.le.3) then

        call copy3d(0,ni+1,0,nj+1,1,nk,up,u)
        call copy3d(0,ni+1,0,nj+1,1,nk,vp,v)
        call copy3d(0,ni+1,0,nj+1,1,nk,wp,w)

        call copy3d(0,ni+1,0,nj+1,1,nk,ppp,pp)

        call copy3d(0,ni+1,0,nj+1,1,nk,ptpp,ptp)

        call copy3d(0,ni+1,0,nj+1,1,nk,qvp,qv)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            call copy4d(0,ni+1,0,nj+1,1,nk,1,nqw,qwtrp,qwtr)

            if(abs(cphopt).eq.4) then

              call copy4d(0,ni+1,0,nj+1,1,nk,1,nnw,nwtrp,nwtr)

            end if

            call setcst4d(0,ni+1,0,nj+1,1,2,1,nqw,0.e0,prwtr)

          end if

          if(abs(cphopt).ge.2) then

            call copy4d(0,ni+1,0,nj+1,1,nk,1,nqi,qicep,qice)

            if(abs(cphopt).eq.2) then

              call copy4d(0,ni+1,0,nj+1,1,nk,1,1,nicep,nice)

            else

              call copy4d(0,ni+1,0,nj+1,1,nk,1,nni,nicep,nice)

            end if

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

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            call setcst4d(0,ni+1,0,nj+1,1,2,1,nqw,0.e0,prwtr)

          end if

          if(abs(cphopt).ge.2) then

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

      end subroutine s_setvar3d

!-----7--------------------------------------------------------------7--

      end module m_setvar3d
