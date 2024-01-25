!***********************************************************************
      module m_smalldrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/02/01
!     Modification: 2004/04/15, 2004/06/10, 2005/01/31, 2006/01/10,
!                   2006/04/03, 2006/06/21, 2006/11/06, 2007/05/07,
!                   2007/07/30, 2008/05/02, 2008/06/09, 2008/08/25,
!                   2008/12/11, 2009/02/27, 2009/03/23, 2011/09/22,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the small time steps
!     integration.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_copy3d
      use m_getiname
      use m_heve
      use m_hevi
      use m_phy2cnt

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: smalldrv, s_smalldrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface smalldrv

        module procedure s_smalldrv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_smalldrv(fpgwmopt,fpimpopt,nsstp,dts,gtinc,          &
     &                      ni,nj,nk,j31,j32,jcb,jcb8u,jcb8v,jcb8w,     &
     &                      mf,mf8u,mf8v,rmf,rmf8u,rmf8v,ubr,vbr,       &
     &                      ptbr,rbr,rst,rst8u,rst8v,rst8w,rcsq,        &
     &                      up,vp,wp,ppp,ptpp,ufrc,vfrc,wfrc,pfrc,ptfrc,&
     &                      ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy,    &
     &                      ptcpx,ptcpy,ugpv,utd,vgpv,vtd,wgpv,wtd,     &
     &                      ppgpv,pptd,ptpgpv,ptptd,uf,vf,wf,ppf,ptpf,  &
     &                      wc,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

      integer, intent(in) :: fpimpopt
                       ! Formal parameter of unique index of impopt

      integer, intent(in) :: nsstp
                       ! Number of small time steps integration

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(in) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

      real, intent(in) :: rcsq(0:ni+1,0:nj+1,1:nk)
                       ! rbr x sound wave speed squared

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(in) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(in) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

      real, intent(in) :: pfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

      real, intent(in) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in potential temperature equation

      real, intent(in) :: ucpx(1:nj,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on west and east boundary

      real, intent(in) :: ucpy(1:ni,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on south and north boundary

      real, intent(in) :: vcpx(1:nj,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on west and east boundary

      real, intent(in) :: vcpy(1:ni,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on south and north boundary

      real, intent(in) :: wcpx(1:nj,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on west and east boundary

      real, intent(in) :: wcpy(1:ni,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on south and north boundary

      real, intent(in) :: pcpx(1:nj,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on west and east boundary

      real, intent(in) :: pcpy(1:ni,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on south and north boundary

      real, intent(in) :: ptcpx(1:nj,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on west and east boundary

      real, intent(in) :: ptcpy(1:ni,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on south and north boundary

      real, intent(in) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(in) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(in) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

      real, intent(in) :: ppgpv(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, intent(in) :: pptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! pressure perturbation of GPV data

      real, intent(in) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(in) :: ptptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

! Output variables

      real, intent(out) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(out) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(out) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

      real, intent(out) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

      real, intent(out) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

! Internal shared variables

      integer gwmopt   ! Option for gravity wave mode integration
      integer impopt   ! Option for vertical implicit method

      integer isstp    ! Index of small time steps integration

      real, intent(inout) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpgwmopt,gwmopt)
      call getiname(fpimpopt,impopt)

! -----

! Copy the past variables to the future variables.

      call copy3d(0,ni+1,0,nj+1,1,nk,up,uf)
      call copy3d(0,ni+1,0,nj+1,1,nk,vp,vf)
      call copy3d(0,ni+1,0,nj+1,1,nk,wp,wf)
      call copy3d(0,ni+1,0,nj+1,1,nk,ppp,ppf)

      if(gwmopt.eq.1) then

        call copy3d(0,ni+1,0,nj+1,1,nk,ptpp,ptpf)

      end if

! -----

! Calculate the zeta components of contravariant velocity.

      call s_phy2cnt(idsthopt,idtrnopt,idmpopt,idmfcopt,idoneopt,       &
     &               ni,nj,nk,j31,j32,jcb8w,mf,uf,vf,wf,wc,             &
     &               tmp1,tmp2,tmp3)

! -----

!! The loop for the small time steps integration.

      do isstp=1,nsstp

! Perform the horizontally explicit and vertically explicit method.

        if(impopt.eq.0) then

          call heve(idgwmopt,idbuyopt,isstp,dts,gtinc,                  &
     &             ni,nj,nk,j31,j32,jcb,jcb8u,jcb8v,jcb8w,              &
     &             mf,mf8u,mf8v,rmf,rmf8u,rmf8v,ubr,vbr,ptbr,rbr,       &
     &             rst,rst8u,rst8v,rst8w,rcsq,ufrc,vfrc,wfrc,pfrc,ptfrc,&
     &             ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy, &
     &             ugpv,utd,vgpv,vtd,wgpv,wtd,ppgpv,pptd,ptpgpv,ptptd,  &
     &             uf,vf,wf,wc,ppf,ptpf,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

! -----

! Perform the horizontally explicit and vertically implicit method.

        else if(impopt.ge.1) then

          call hevi(idgwmopt,idbuyopt,isstp,dts,gtinc,                  &
     &             ni,nj,nk,j31,j32,jcb,jcb8u,jcb8v,jcb8w,              &
     &             mf,mf8u,mf8v,rmf,rmf8u,rmf8v,ubr,vbr,ptbr,rbr,       &
     &             rst,rst8u,rst8v,rst8w,rcsq,ufrc,vfrc,wfrc,pfrc,ptfrc,&
     &             ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy, &
     &             ugpv,utd,vgpv,vtd,wgpv,wtd,ppgpv,pptd,ptpgpv,ptptd,  &
     &             uf,vf,wf,wc,ppf,ptpf,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

        end if

! -----

      end do

!! -----

      end subroutine s_smalldrv

!-----7--------------------------------------------------------------7--

      end module m_smalldrv
