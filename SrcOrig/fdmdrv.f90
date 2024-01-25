!***********************************************************************
      module m_fdmdrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/04/03
!     Modification: 2006/06/21, 2006/07/21, 2006/11/06, 2007/01/20,
!                   2007/05/07, 2007/07/30, 2008/05/02, 2008/08/25,
!                   2008/12/11, 2009/02/27, 2009/03/23, 2010/04/01,
!                   2011/08/18, 2011/09/22, 2011/11/10

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the finite difference time
!     steps integration.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_comkind
      use m_forcedrv
      use m_smalldrv
      use m_steps

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: fdmdrv, s_fdmdrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface fdmdrv

        module procedure s_fdmdrv

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
      subroutine s_fdmdrv(fmois,ksp0,nsstp,ctime,                       &
     &                dtb,dts,nggdmp,ngrdmp,gtinc,atinc,rtinc,          &
     &                ni,nj,nk,nqw,nnw,nqi,nni,nqa,zph,                 &
     &                j31,j32,jcb,jcb8u,jcb8v,jcb8w,mf,mf8u,mf8v,       &
     &                rmf,rmf8u,rmf8v,fc,ubr,vbr,pbr,ptbr,qvbr,rbr,     &
     &                rst,rst8u,rst8v,rst8w,rcsq,u,up,v,vp,w,wp,        &
     &                pp,ppp,ptp,ptpp,qv,qvp,qwtr,qwtrp,nwtr,nwtrp,     &
     &                qice,qicep,nice,nicep,qcwtr,qcwtrp,qcice,qcicep,  &
     &                qasl,qaslp,qt,qtp,ucpx,ucpy,vcpx,vcpy,            &
     &                wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,      &
     &                qwcpx,qwcpy,nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,  &
     &                qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,          &
     &                qtcpx,qtcpy,tkecpx,tkecpy,rbcx,rbcy,rbcxy,rbct,   &
     &                ugpv,utd,vgpv,vtd,wgpv,wtd,ppgpv,pptd,            &
     &                ptpgpv,ptptd,qvgpv,qvtd,qwgpv,qwtd,qigpv,qitd,    &
     &                qagpv,qatd,urdr,vrdr,wrdr,qwrdr,qwrtd,qirdr,qirtd,&
     &                qall,qallp,tke,tkep,ufrc,vfrc,pfrc,ptfrc,qvfrc,   &
     &                uf,vf,wf,wc,ppf,ptpf,qvf,qwtrf,nwtrf,qicef,nicef, &
     &                qcwtrf,qcicef,qaslf,qtf,tkef,wfrc,tmp1,tmp2,tmp3, &
     &                tmp4,tmp5,tmp6,tmp7)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      integer, intent(in) :: nsstp
                       ! Number of small time steps integration

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

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

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient for GPV

      real, intent(in) :: ngrdmp(1:2)
                       ! Analysis nudging damping coefficient for radar

      real, intent(in) :: gtinc
                       ! Lapse of forecast time
                       ! from GPV data reading

      real, intent(in) :: atinc
                       ! Lapse of forecast time
                       ! from aerosol data reading

      real, intent(in) :: rtinc(1:2)
                       ! Lapse of forecast time
                       ! from radar data reading

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

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

      real, intent(in) :: fc(0:ni+1,0:nj+1,1:2)
                       ! 0.25 x Coriolis parameters

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

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

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at present

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at present

      real, intent(in) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(in) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at present

      real, intent(in) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

      real, intent(in) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(in) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at present

      real, intent(in) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(in) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at present

      real, intent(in) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(in) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at present

      real, intent(in) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(in) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(in) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(in) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at present

      real, intent(in) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

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

      real, intent(in) :: qvcpx(1:nj,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on west and east boundary

      real, intent(in) :: qvcpy(1:ni,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on south and north boundary

      real, intent(in) :: qwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on west and east boundary

      real, intent(in) :: qwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on south and north boundary

      real, intent(in) :: nwcpx(1:nj,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on west and east boundary

      real, intent(in) :: nwcpy(1:ni,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on south and north boundary

      real, intent(in) :: qicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on west and east boundary

      real, intent(in) :: qicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on south and north boundary

      real, intent(in) :: nicpx(1:nj,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on west and east boundary

      real, intent(in) :: nicpy(1:ni,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on south and north boundary

      real, intent(in) :: qcwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on west and east boundary

      real, intent(in) :: qcwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on south and north boundary

      real, intent(in) :: qcicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on west and east boundary

      real, intent(in) :: qcicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on south and north boundary

      real, intent(in) :: qacpx(1:nj,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on west and east boundary

      real, intent(in) :: qacpy(1:ni,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on south and north boundary

      real, intent(in) :: qtcpx(1:nj,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on west and east boundary

      real, intent(in) :: qtcpy(1:ni,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on south and north boundary

      real, intent(in) :: tkecpx(1:nj,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on west and east boundary

      real, intent(in) :: tkecpy(1:ni,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on south and north boundary

      real, intent(in) :: rbcx(1:ni)
                       ! Relaxed lateral sponge damping coefficients
                       ! in x direction

      real, intent(in) :: rbcy(1:nj)
                       ! Relaxed lateral sponge damping coefficients
                       ! in y direction

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

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

      real, intent(in) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(in) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

      real, intent(in) :: qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of GPV data at marked time

      real, intent(in) :: qwtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of water hydrometeor of GPV data

      real, intent(in) :: qigpv(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of GPV data at marked time

      real, intent(in) :: qitd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of ice hydrometeor of GPV data

      real, intent(in) :: qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Mixing ratio of aerosol data at marked time

      real, intent(in) :: qatd(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Time tendency of mixing ratio of aerosol data

      real, intent(in) :: urdr(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of radar data
                       ! at marked time

      real, intent(in) :: vrdr(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of radar data
                       ! at marked time

      real, intent(in) :: wrdr(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of radar data
                       ! at marked time

      real, intent(in) :: qwrdr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of radar data
                       ! at marked time

      real, intent(in) :: qwrtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of
                       ! water hydrometeor of radar data

      real, intent(in) :: qirdr(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of radar data
                       ! at marked time

      real, intent(in) :: qirtd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of
                       ! ice hydrometeor of radar data

      real, intent(in) :: qall(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at present

      real, intent(in) :: qallp(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at past

! Input and output variables

      real, intent(inout) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(inout) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: pfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

      real, intent(inout) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in potential temperature equation

      real, intent(inout) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

! Output variables

      real, intent(out) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(out) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(out) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

      real, intent(out) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

      real, intent(out) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

      real, intent(out) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(out) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(out) :: qwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at future

      real, intent(out) :: nwtrf(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at future

      real, intent(out) :: qicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at future

      real, intent(out) :: nicef(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at future

      real, intent(out) :: qcwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at future

      real, intent(out) :: qcicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at future

      real, intent(out) :: qaslf(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at future

      real, intent(out) :: qtf(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at future

      real, intent(out) :: tkef(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at future

! Internal shared variables

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

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

      real, intent(inout) :: tmp7(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Remark

!     uf,vf,wf,ppf,ptpf,qvf: These variables are also temporary.

!-----7--------------------------------------------------------------7--

! Perform the forcing terms calculation in the large time steps.

      call forcedrv(idtubopt,fmois,ksp0,                                &
     &              ctime,dtb,nggdmp,ngrdmp,gtinc,atinc,rtinc,          &
     &              ni,nj,nk,nqw,nnw,nqi,nni,nqa,zph,j31,j32,jcb,       &
     &              jcb8u,jcb8v,jcb8w,mf,mf8u,mf8v,rmf,rmf8u,rmf8v,fc,  &
     &              ubr,vbr,pbr,ptbr,qvbr,rbr,rst,rst8u,rst8v,rst8w,    &
     &              u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,qv,qvp,              &
     &              qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,        &
     &              qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,        &
     &              rbcx,rbcy,rbcxy,rbct,ugpv,utd,vgpv,vtd,wgpv,wtd,    &
     &              ppgpv,pptd,ptpgpv,ptptd,qvgpv,qvtd,                 &
     &              qwgpv,qwtd,qigpv,qitd,qagpv,qatd,urdr,vrdr,wrdr,    &
     &              qwrdr,qwrtd,qirdr,qirtd,qall,qallp,                 &
     &              tke,tkep,ufrc,vfrc,pfrc,ptfrc,qvfrc,wfrc,           &
     &              qwtrf,nwtrf,qicef,nicef,qcwtrf,qcicef,              &
     &              qaslf,qtf,tkef,wc,uf,vf,wf,ppf,ptpf,qvf,            &
     &              tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)

! -----

! Perform the small time steps integration.

      call smalldrv(idgwmopt,idimpopt,nsstp,dts,gtinc,ni,nj,nk,j31,j32, &
     &              jcb,jcb8u,jcb8v,jcb8w,mf,mf8u,mf8v,rmf,rmf8u,rmf8v, &
     &              ubr,vbr,ptbr,rbr,rst,rst8u,rst8v,rst8w,rcsq,        &
     &              up,vp,wp,ppp,ptpp,ufrc,vfrc,wfrc,pfrc,ptfrc,        &
     &              ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy,ptcpx,ptcpy,&
     &              ugpv,utd,vgpv,vtd,wgpv,wtd,ppgpv,pptd,ptpgpv,ptptd, &
     &              uf,vf,wf,ppf,ptpf,wc,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

! -----

! Solve the scalar variables to the next time step.

      call steps(idgpvvar,idexbvar,idexbopt,idgwmopt,idadvopt,          &
     &           idcphopt,idhaiopt,idqcgopt,idaslopt,idtrkopt,idtubopt, &
     &           fmois,dtb,gtinc,atinc,ni,nj,nk,nqw,nnw,nqi,nni,nqa,    &
     &           qvbr,rst,ptp,ptpp,qv,qvp,qwtr,qwtrp,nwtr,nwtrp,        &
     &           qice,qicep,nice,nicep,qcwtr,qcwtrp,qcice,qcicep,       &
     &           qasl,qaslp,qt,qtp,tke,tkep,ptcpx,ptcpy,                &
     &           qvcpx,qvcpy,qwcpx,qwcpy,nwcpx,nwcpy,qicpx,qicpy,       &
     &           nicpx,nicpy,qcwcpx,qcwcpy,qcicpx,qcicpy,               &
     &           qacpx,qacpy,qtcpx,qtcpy,tkecpx,tkecpy,                 &
     &           ptpgpv,ptptd,qvgpv,qvtd,qwgpv,qwtd,qigpv,qitd,         &
     &           qagpv,qatd,ptfrc,qvfrc,qwtrf,nwtrf,qicef,nicef,        &
     &           qcwtrf,qcicef,qaslf,qtf,tkef,ptpf,qvf,tmp1)

! -----

      end subroutine s_fdmdrv

!-----7--------------------------------------------------------------7--

      end module m_fdmdrv
