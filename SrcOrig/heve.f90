!***********************************************************************
      module m_heve
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/12/17
!     Modification: 1999/12/20, 2000/01/17, 2000/03/17, 2000/12/18,
!                   2001/01/15, 2001/02/24, 2001/03/13, 2001/05/29,
!                   2001/06/06, 2001/06/29, 2001/07/13, 2001/08/07,
!                   2001/09/13, 2001/11/20, 2002/04/02, 2002/07/23,
!                   2002/08/15, 2003/01/04, 2003/03/21, 2003/04/30,
!                   2003/05/19, 2003/11/05, 2003/11/28, 2003/12/12,
!                   2003/12/26, 2004/04/15, 2004/05/31, 2004/06/10,
!                   2004/07/01, 2005/01/31, 2006/04/03, 2006/06/21,
!                   2006/11/06, 2007/05/07, 2008/05/02, 2008/08/25,
!                   2008/12/11, 2009/01/30, 2009/02/27, 2009/03/23,
!                   2011/09/22, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the horizontally explicit and
!     vertically explicit method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advbspe
      use m_advbspt
      use m_buoywse
      use m_comindx
      use m_diverpe
      use m_getiname
      use m_pgrad
      use m_steppe
      use m_steppts
      use m_stepuv
      use m_stepwe

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: heve, s_heve

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface heve

        module procedure s_heve

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
      subroutine s_heve(fpgwmopt,fpbuyopt,isstp,dts,gtinc,              &
     &                  ni,nj,nk,j31,j32,jcb,jcb8u,jcb8v,jcb8w,         &
     &                  mf,mf8u,mf8v,rmf,rmf8u,rmf8v,ubr,vbr,ptbr,rbr,  &
     &                  rst,rst8u,rst8v,rst8w,rcsq,ufrc,vfrc,wfrc,      &
     &                  pfrc,ptfrc,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,       &
     &                  pcpx,pcpy,ptcpx,ptcpy,ugpv,utd,vgpv,vtd,        &
     &                  wgpv,wtd,ppgpv,pptd,ptpgpv,ptptd,uf,vf,wf,wc,   &
     &                  ppf,ptpf,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

      integer, intent(in) :: fpbuyopt
                       ! Formal parameter of unique index of buyopt

      integer, intent(in) :: isstp
                       ! Index of small time steps integration

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

! Input and output variables

      real, intent(inout) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(inout) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(inout) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

      real, intent(inout) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity
                       ! at future

      real, intent(inout) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

! Internal shared variables

      integer gwmopt   ! Option for gravity wave mode integration
      integer buyopt   ! Option for buoyancy calculation

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
      call getiname(fpbuyopt,buyopt)

! -----

! Calculate the pressure gradient force.

      call pgrad(idtrnopt,idmpopt,idmfcopt,iddivopt,iddx,iddy,iddz,     &
     &           iddxiv,iddyiv,iddziv,dts,ni,nj,nk,j31,j32,jcb,         &
     &           mf,mf8u,mf8v,rmf,rmf8u,rmf8v,rst8u,rst8v,rst8w,        &
     &           uf,vf,wc,ppf,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

! -----

! Solve the x and y components of velocity to the next time step.

      call stepuv(idexbvar,idexbopt,isstp,dts,gtinc,                    &
     &            ni,nj,nk,ubr,vbr,rst8u,rst8v,ufrc,vfrc,               &
     &            tmp1,tmp2,ucpx,ucpy,vcpx,vcpy,ugpv,utd,vgpv,vtd,uf,vf)

! -----

! Calculate the buoyancy in the small time steps.

      if(buyopt.eq.1) then

        call buoywse(idgwmopt,ni,nj,nk,ptbr,rst,rcsq,ppf,ptpf,tmp3,tmp4)

      end if

! -----

!! Solve the potential temperature perturbation.

      if(gwmopt.eq.1) then

! Calculate the base state advection in the potential temperature
! eqation.

        call advbspt(idgwmopt,iddziv,ni,nj,nk,ptbr,rbr,wf,tmp6,tmp1)

! -----

! Solve the potential temperature perturbation to the next time step.

        call steppts(idexbvar,idexbopt,isstp,dts,gtinc,ni,nj,nk,        &
     &               rst,ptfrc,tmp6,ptcpx,ptcpy,ptpgpv,ptptd,ptpf)

! -----

      end if

!! -----

! Solve the z components of velocity to the next time step.

      call stepwe(idexbvar,idexbopt,isstp,dts,gtinc,ni,nj,nk,j31,j32,   &
     &            jcb8w,mf,rst8w,uf,vf,wfrc,tmp3,wcpx,wcpy,wgpv,wtd,    &
     &            wf,wc,tmp4,tmp5,tmp6)

! -----

! Calculate the base state pressure advection and the divergence in the
! pressure equation.

      call diverpe(ni,nj,nk,jcb8u,jcb8v,jcb8w,mf,rmf,rmf8u,rmf8v,rcsq,  &
     &             uf,vf,wc,tmp2,tmp3,tmp4,tmp5)

      call advbspe(ni,nj,nk,rst,wf,tmp2,tmp1)

! -----

! Solve the pressure perturbation to the next time step.

      call steppe(idexbvar,idexbopt,isstp,dts,gtinc,ni,nj,nk,jcb,pfrc,  &
     &            tmp1,pcpx,pcpy,ppgpv,pptd,ppf)

! -----

      end subroutine s_heve

!-----7--------------------------------------------------------------7--

      end module m_heve
