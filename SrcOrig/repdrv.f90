!***********************************************************************
      module m_repdrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/07/01
!     Modification: 2005/08/05, 2005/10/05, 2006/01/10, 2006/02/13,
!                   2006/04/03, 2006/07/21, 2007/01/20, 2007/05/07,
!                   2007/07/30, 2007/11/26, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/01/30, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2011/11/10, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the repositon.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_getiname
      use m_inichar
      use m_repsit3d
      use m_repsit2d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: repdrv, s_repdrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface repdrv

        module procedure s_repdrv

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
      subroutine s_repdrv(fpdmpvar,                                     &
     &                  fpwbc,fpebc,fpsbc,fpnbc,fpsfcopt,fpadvopt,      &
     &                  fpcphopt,fphaiopt,fpqcgopt,fpaslopt,fptrkopt,   &
     &                  fptubopt,fpdiaopt,fmois,istr,iend,jstr,jend,    &
     &                  di,dj,istrb,iendb,jstrb,jendb,dib,djb,          &
     &                  ni,nj,nk,nqw,nnw,nqi,nni,nqa,nund,              &
     &                  ni_rst,nj_rst,ubr,vbr,pbr,ptbr,qvbr,            &
     &                  u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,qv,qvp,          &
     &                  qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,    &
     &                  qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,    &
     &                  tke,tkep,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,         &
     &                  pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,              &
     &                  qwcpx,qwcpy,nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,&
     &                  qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,        &
     &                  qtcpx,qtcpy,tkecpx,tkecpy,maxvl,                &
     &                  prwtr,price,pdia,z0m,z0h,tund,tundp,            &
     &                  ubr_rst,vbr_rst,pbr_rst,ptbr_rst,qvbr_rst,      &
     &                  u_rst,up_rst,v_rst,vp_rst,w_rst,wp_rst,         &
     &                  pp_rst,ppp_rst,ptp_rst,ptpp_rst,qv_rst,qvp_rst, &
     &                  qwtr_rst,qwtrp_rst,nwtr_rst,nwtrp_rst,          &
     &                  qice_rst,qicep_rst,nice_rst,nicep_rst,          &
     &                  qcwtr_rst,qcwtrp_rst,qcice_rst,qcicep_rst,      &
     &                  qasl_rst,qaslp_rst,qt_rst,qtp_rst,              &
     &                  tke_rst,tkep_rst,ucpx_rst,ucpy_rst,             &
     &                  vcpx_rst,vcpy_rst,wcpx_rst,wcpy_rst,            &
     &                  pcpx_rst,pcpy_rst,ptcpx_rst,ptcpy_rst,          &
     &                  qvcpx_rst,qvcpy_rst,qwcpx_rst,qwcpy_rst,        &
     &                  nwcpx_rst,nwcpy_rst,qicpx_rst,qicpy_rst,        &
     &                  nicpx_rst,nicpy_rst,qcwcpx_rst,qcwcpy_rst,      &
     &                  qcicpx_rst,qcicpy_rst,qacpx_rst,qacpy_rst,      &
     &                  qtcpx_rst,qtcpy_rst,tkecpx_rst,tkecpy_rst,      &
     &                  maxvl_rst,prwtr_rst,price_rst,pdia_rst,         &
     &                  z0m_rst,z0h_rst,tund_rst,tundp_rst)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

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

      integer, intent(in) :: fpdiaopt
                       ! Formal parameter of unique index of diaopt

      integer, intent(in) :: istr
                       ! Minimum do loops index in x direction

      integer, intent(in) :: iend
                       ! Maximum do loops index in x direction

      integer, intent(in) :: jstr
                       ! Minimum do loops index in y direction

      integer, intent(in) :: jend
                       ! Maximum do loops index in y direction

      integer, intent(in) :: di
                       ! Differential index to istr

      integer, intent(in) :: dj
                       ! Differential index to jstr

      integer, intent(in) :: istrb
                       ! Minimum do loops index in x direction
                       ! of lateral boundary

      integer, intent(in) :: iendb
                       ! Maximum do loops index in x direction
                       ! of lateral boundary

      integer, intent(in) :: jstrb
                       ! Minimum do loops index in y direction
                       ! of lateral boundary

      integer, intent(in) :: jendb
                       ! Maximum do loops index in y direction
                       ! of lateral boundary

      integer, intent(in) :: dib
                       ! Differential index to istrb

      integer, intent(in) :: djb
                       ! Differential index to jstrb

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

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: ni_rst
                       ! Restructed files dimension in x direction

      integer, intent(in) :: nj_rst
                       ! Restructed files dimension in y direction

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

      real, intent(in) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(in) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

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

      real, intent(in) :: maxvl(0:ni+1,0:nj+1,1:nk)
                       ! Maximum instantaneous wind velocity

      real, intent(in) :: prwtr(0:ni+1,0:nj+1,1:2,1:nqw)
                       ! Precipitation and accumulation for water

      real, intent(in) :: price(0:ni+1,0:nj+1,1:2,1:nqi)
                       ! Precipitation and accumulation for ice

      real, intent(in) :: pdia(0:ni+1,0:nj+1,1:nk)
                       ! Diabatic value

      real, intent(in) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(in) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(in) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at present

      real, intent(in) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

! Input and output variables

      real, intent(inout) :: ubr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ubr in restructed domain

      real, intent(inout) :: vbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! vbr in restructed domain

      real, intent(inout) :: pbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! pbr in restructed domain

      real, intent(inout) :: ptbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ptbr in restructed domain

      real, intent(inout) :: qvbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qvbr in restructed domain

      real, intent(inout) :: u_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! u in restructed domain

      real, intent(inout) :: up_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! up in restructed domain

      real, intent(inout) :: v_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! v in restructed domain

      real, intent(inout) :: vp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! vp in restructed domain

      real, intent(inout) :: w_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! w in restructed domain

      real, intent(inout) :: wp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! wp in restructed domain

      real, intent(inout) :: pp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! pp in restructed domain

      real, intent(inout) :: ppp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ppp in restructed domain

      real, intent(inout) :: ptp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ptp in restructed domain

      real, intent(inout) :: ptpp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! ptpp in restructed domain

      real, intent(inout) :: qv_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qv in restructed domain

      real, intent(inout) :: qvp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qvp in restructed domain

      real, intent(inout) :: qwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qwtr in restructed domain

      real, intent(inout) :: qwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qwtrp in restructed domain

      real, intent(inout) :: nwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nnw)
                       ! nwce in restructed domain

      real, intent(inout) :: nwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nnw)
                       ! nwcep in restructed domain

      real, intent(inout) :: qice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qice in restructed domain

      real, intent(inout) :: qicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qicep in restructed domain

      real, intent(inout) :: nice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nni)
                       ! nice in restructed domain

      real, intent(inout) :: nicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nni)
                       ! nicep in restructed domain

      real, intent(inout) ::                                            &
     &                   qcwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qcwtr in restructed domain

      real, intent(inout) ::                                            &
     &                   qcwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw)
                       ! qcwtrp in restructed domain

      real, intent(inout) ::                                            &
     &                   qcice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qcice in restructed domain

      real, intent(inout) ::                                            &
     &                   qcicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi)
                       ! qcicep in restructed domain

      real, intent(inout) ::                                            &
     &                   qasl_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqa(0))
                       ! qasl in restructed domain

      real, intent(inout) ::                                            &
     &                   qaslp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqa(0))
                       ! qaslp in restructed domain

      real, intent(inout) :: qt_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qt in restructed domain

      real, intent(inout) :: qtp_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! qtp in restructed domain

      real, intent(inout) :: tke_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! tke in restructed domain

      real, intent(inout) :: tkep_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! tkep in restructed domain

      real, intent(inout) :: ucpx_rst(1:nj_rst,1:nk,1:2)
                       ! ucpx in restructed domain

      real, intent(inout) :: ucpy_rst(1:ni_rst,1:nk,1:2)
                       ! ucpy in restructed domain

      real, intent(inout) :: vcpx_rst(1:nj_rst,1:nk,1:2)
                       ! vcpx in restructed domain

      real, intent(inout) :: vcpy_rst(1:ni_rst,1:nk,1:2)
                       ! vcpy in restructed domain

      real, intent(inout) :: wcpx_rst(1:nj_rst,1:nk,1:2)
                       ! wcpx in restructed domain

      real, intent(inout) :: wcpy_rst(1:ni_rst,1:nk,1:2)
                       ! wcpy in restructed domain

      real, intent(inout) :: pcpx_rst(1:nj_rst,1:nk,1:2)
                       ! pcpx in restructed domain

      real, intent(inout) :: pcpy_rst(1:ni_rst,1:nk,1:2)
                       ! pcpy in restructed domain

      real, intent(inout) :: ptcpx_rst(1:nj_rst,1:nk,1:2)
                       ! ptcpx in restructed domain

      real, intent(inout) :: ptcpy_rst(1:ni_rst,1:nk,1:2)
                       ! ptcpy in restructed domain

      real, intent(inout) :: qvcpx_rst(1:nj_rst,1:nk,1:2)
                       ! qvcpx in restructed domain

      real, intent(inout) :: qvcpy_rst(1:ni_rst,1:nk,1:2)
                       ! qvcpy in restructed domain

      real, intent(inout) :: qwcpx_rst(1:nj_rst,1:nk,1:2,1:nqw)
                       ! qwcpx in restructed domain

      real, intent(inout) :: qwcpy_rst(1:ni_rst,1:nk,1:2,1:nqw)
                       ! qwcpy in restructed domain

      real, intent(inout) :: nwcpx_rst(1:nj_rst,1:nk,1:2,1:nnw)
                       ! nwcpx in restructed domain

      real, intent(inout) :: nwcpy_rst(1:ni_rst,1:nk,1:2,1:nnw)
                       ! nwcpy in restructed domain

      real, intent(inout) :: qicpx_rst(1:nj_rst,1:nk,1:2,1:nqi)
                       ! qicpx in restructed domain

      real, intent(inout) :: qicpy_rst(1:ni_rst,1:nk,1:2,1:nqi)
                       ! qicpy in restructed domain

      real, intent(inout) :: nicpx_rst(1:nj_rst,1:nk,1:2,1:nni)
                       ! nicpx in restructed domain

      real, intent(inout) :: nicpy_rst(1:ni_rst,1:nk,1:2,1:nni)
                       ! nicpy in restructed domain

      real, intent(inout) :: qcwcpx_rst(1:nj_rst,1:nk,1:2,1:nqw)
                       ! qcwcpx in restructed domain

      real, intent(inout) :: qcwcpy_rst(1:ni_rst,1:nk,1:2,1:nqw)
                       ! qcwcpy in restructed domain

      real, intent(inout) :: qcicpx_rst(1:nj_rst,1:nk,1:2,1:nqi)
                       ! qcicpx in restructed domain

      real, intent(inout) :: qcicpy_rst(1:ni_rst,1:nk,1:2,1:nqi)
                       ! qcicpy in restructed domain

      real, intent(inout) :: qacpx_rst(1:nj_rst,1:nk,1:2,1:nqa(0))
                       ! qacpx in restructed domain

      real, intent(inout) :: qacpy_rst(1:ni_rst,1:nk,1:2,1:nqa(0))
                       ! qacpy in restructed domain

      real, intent(inout) :: qtcpx_rst(1:nj_rst,1:nk,1:2)
                       ! qtcpx in restructed domain

      real, intent(inout) :: qtcpy_rst(1:ni_rst,1:nk,1:2)
                       ! qtcpy in restructed domain

      real, intent(inout) :: tkecpx_rst(1:nj_rst,1:nk,1:2)
                       ! tkecpx in restructed domain

      real, intent(inout) :: tkecpy_rst(1:ni_rst,1:nk,1:2)
                       ! tkecpy in restructed domain

      real, intent(inout) :: maxvl_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! maxvl in restructed domain

      real, intent(inout) :: prwtr_rst(0:ni_rst+1,0:nj_rst+1,1:2,1:nqw)
                       ! prwtr in restructed domain

      real, intent(inout) :: price_rst(0:ni_rst+1,0:nj_rst+1,1:2,1:nqi)
                       ! price in restructed domain

      real, intent(inout) :: pdia_rst(0:ni_rst+1,0:nj_rst+1,1:nk)
                       ! pdia in restructed domain

      real, intent(inout) :: z0m_rst(0:ni_rst+1,0:nj_rst+1)
                       ! z0m in restructed domain

      real, intent(inout) :: z0h_rst(0:ni_rst+1,0:nj_rst+1)
                       ! z0h in restructed domain

      real, intent(inout) :: tund_rst(0:ni_rst+1,0:nj_rst+1,1:nund)
                       ! tund in restructed domain

      real, intent(inout) :: tundp_rst(0:ni_rst+1,0:nj_rst+1,1:nund)
                       ! tundp in restructed domain

! Internal shared variables

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer sfcopt   ! Optoin for surface physics
      integer advopt   ! Optoin for advection scheme
      integer cphopt   ! Optoin for cloud micro physics
      integer haiopt   ! Optoin for additional hail processes
      integer qcgopt   ! Optoin for charging distribution
      integer aslopt   ! Optoin for aerosol processes
      integer trkopt   ! Optoin for mixing ratio tracking
      integer tubopt   ! Optoin for turbulent mixing
      integer diaopt   ! Optoin for diabatic calculation

      integer n        ! Array index in 4th direction

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(dmpvar)

! -----

! Get the required namelist variables.

      call getcname(fpdmpvar,dmpvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpsfcopt,sfcopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)
      call getiname(fpdiaopt,diaopt)

! -----

!!!!! Reposition the restructed variables form original restart
!!!!! variables.

      if(advopt.le.3) then

!!!! Reposition the restructed variables in the case the centered
!!!! advection scheme is performed.

! For base state variables.

        call repsit3d('ox',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,ubr,ubr_rst)

        call repsit3d('xo',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,vbr,vbr_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,pbr,pbr_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,ptbr,ptbr_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,qvbr,qvbr_rst)

! -----

! For the velocity, the pressure and the potential temperature.

        call repsit3d('ox',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,u,u_rst)

        call repsit3d('ox',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,up,up_rst)

        call repsit3d('xo',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,v,v_rst)

        call repsit3d('xo',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,vp,vp_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,w,w_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,wp,wp_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,pp,pp_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,ppp,ppp_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,ptp,ptp_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,ptpp,ptpp_rst)

! -----

! For the boundary variables.

        if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

          call repsit2d('ox',istrb,iendb,jstrb,jendb,dib,djb,ni,nj,nk,  &
     &                  ni_rst,nj_rst,ucpx,ucpy,ucpx_rst,ucpy_rst)

          call repsit2d('xo',istrb,iendb,jstrb,jendb,dib,djb,ni,nj,nk,  &
     &                  ni_rst,nj_rst,vcpx,vcpy,vcpx_rst,vcpy_rst)

          call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,ni,nj,nk,  &
     &                  ni_rst,nj_rst,wcpx,wcpy,wcpx_rst,wcpy_rst)

          call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,ni,nj,nk,  &
     &                  ni_rst,nj_rst,pcpx,pcpy,pcpx_rst,pcpy_rst)

          call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,ni,nj,nk,  &
     &                  ni_rst,nj_rst,ptcpx,ptcpy,ptcpx_rst,ptcpy_rst)

        end if

! -----

!!! For the water substance.

        if(fmois(1:5).eq.'moist') then

! For the water vapor mixing ratio.

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,        &
     &                  ni_rst,nj_rst,qv,qv_rst)

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,        &
     &                  ni_rst,nj_rst,qvp,qvp_rst)

          if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

            call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,         &
     &                    ni,nj,nk,ni_rst,nj_rst,qvcpx,qvcpy,           &
     &                    qvcpx_rst,qvcpy_rst)

          end if

! -----

!! For the water and ice hydrometeor.

! For the bulk categories.

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qwtr(0,0,1,1),qwtr_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qwtr(0,0,1,2),qwtr_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qwtrp(0,0,1,1),qwtrp_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qwtrp(0,0,1,2),qwtrp_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        prwtr(0,0,1,1),prwtr_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        prwtr(0,0,1,2),prwtr_rst(0,0,1,2))

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qwcpx(1,1,1,1),qwcpy(1,1,1,1),          &
     &                          qwcpx_rst(1,1,1,1),qwcpy_rst(1,1,1,1))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qwcpx(1,1,1,2),qwcpy(1,1,1,2),          &
     &                          qwcpx_rst(1,1,1,2),qwcpy_rst(1,1,1,2))

              end if

            end if

            if(abs(cphopt).eq.4) then

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nwtr(0,0,1,1),nwtr_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nwtr(0,0,1,2),nwtr_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nwtrp(0,0,1,1),nwtrp_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nwtrp(0,0,1,2),nwtrp_rst(0,0,1,2))

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nwcpx(1,1,1,1),nwcpy(1,1,1,1),          &
     &                          nwcpx_rst(1,1,1,1),nwcpy_rst(1,1,1,1))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nwcpx(1,1,1,2),nwcpy(1,1,1,2),          &
     &                          nwcpx_rst(1,1,1,2),nwcpy_rst(1,1,1,2))

              end if

            end if

            if(abs(cphopt).ge.2) then

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qice(0,0,1,1),qice_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qice(0,0,1,2),qice_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qice(0,0,1,3),qice_rst(0,0,1,3))

              if(haiopt.eq.1) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qice(0,0,1,4),qice_rst(0,0,1,4))

              end if

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qicep(0,0,1,1),qicep_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qicep(0,0,1,2),qicep_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qicep(0,0,1,3),qicep_rst(0,0,1,3))

              if(haiopt.eq.1) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qicep(0,0,1,4),qicep_rst(0,0,1,4))

              end if

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        price(0,0,1,1),price_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        price(0,0,1,2),price_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        price(0,0,1,3),price_rst(0,0,1,3))

              if(haiopt.eq.1) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,2,ni_rst,nj_rst,            &
     &                          price(0,0,1,4),price_rst(0,0,1,4))

              end if

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qicpx(1,1,1,1),qicpy(1,1,1,1),          &
     &                          qicpx_rst(1,1,1,1),qicpy_rst(1,1,1,1))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qicpx(1,1,1,2),qicpy(1,1,1,2),          &
     &                          qicpx_rst(1,1,1,2),qicpy_rst(1,1,1,2))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qicpx(1,1,1,3),qicpy(1,1,1,3),          &
     &                          qicpx_rst(1,1,1,3),qicpy_rst(1,1,1,3))

                if(haiopt.eq.1) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            qicpx(1,1,1,4),qicpy(1,1,1,4),        &
     &                            qicpx_rst(1,1,1,4),qicpy_rst(1,1,1,4))

                end if

              end if

            end if

            if(abs(cphopt).eq.2) then

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nice(0,0,1,1),nice_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nicep(0,0,1,1),nicep_rst(0,0,1,1))

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nicpx(1,1,1,1),nicpy(1,1,1,1),          &
     &                          nicpx_rst(1,1,1,1),nicpy_rst(1,1,1,1))

              end if

            else if(abs(cphopt).ge.3) then

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nice(0,0,1,1),nice_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nice(0,0,1,2),nice_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nice(0,0,1,3),nice_rst(0,0,1,3))

              if(haiopt.eq.1) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          nice(0,0,1,4),nice_rst(0,0,1,4))

              end if

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nicep(0,0,1,1),nicep_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nicep(0,0,1,2),nicep_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nicep(0,0,1,3),nicep_rst(0,0,1,3))

              if(haiopt.eq.1) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          nicep(0,0,1,4),nicep_rst(0,0,1,4))

              end if

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nicpx(1,1,1,1),nicpy(1,1,1,1),          &
     &                          nicpx_rst(1,1,1,1),nicpy_rst(1,1,1,1))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nicpx(1,1,1,2),nicpy(1,1,1,2),          &
     &                          nicpx_rst(1,1,1,2),nicpy_rst(1,1,1,2))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nicpx(1,1,1,3),nicpy(1,1,1,3),          &
     &                          nicpx_rst(1,1,1,3),nicpy_rst(1,1,1,3))

                if(haiopt.eq.1) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            nicpx(1,1,1,4),nicpy(1,1,1,4),        &
     &                            nicpx_rst(1,1,1,4),nicpy_rst(1,1,1,4))

                end if

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qcwtr(0,0,1,1),qcwtr_rst(0,0,1,1))

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qcwtr(0,0,1,2),qcwtr_rst(0,0,1,2))

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qcwtrp(0,0,1,1),qcwtrp_rst(0,0,1,1))

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qcwtrp(0,0,1,2),qcwtrp_rst(0,0,1,2))

              end if

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qcice(0,0,1,1),qcice_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qcice(0,0,1,2),qcice_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qcice(0,0,1,3),qcice_rst(0,0,1,3))

              if(haiopt.eq.1) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qcice(0,0,1,4),qcice_rst(0,0,1,4))

              end if

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qcicep(0,0,1,1),qcicep_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qcicep(0,0,1,2),qcicep_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qcicep(0,0,1,3),qcicep_rst(0,0,1,3))

              if(haiopt.eq.1) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qcicep(0,0,1,4),qcicep_rst(0,0,1,4))

              end if

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                if(qcgopt.eq.2) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcwcpx(1,1,1,1),qcwcpy(1,1,1,1),        &
     &                          qcwcpx_rst(1,1,1,1),qcwcpy_rst(1,1,1,1))

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcwcpx(1,1,1,2),qcwcpy(1,1,1,2),        &
     &                          qcwcpx_rst(1,1,1,2),qcwcpy_rst(1,1,1,2))

                end if

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcicpx(1,1,1,1),qcicpy(1,1,1,1),        &
     &                          qcicpx_rst(1,1,1,1),qcicpy_rst(1,1,1,1))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcicpx(1,1,1,2),qcicpy(1,1,1,2),        &
     &                          qcicpx_rst(1,1,1,2),qcicpy_rst(1,1,1,2))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcicpx(1,1,1,3),qcicpy(1,1,1,3),        &
     &                          qcicpx_rst(1,1,1,3),qcicpy_rst(1,1,1,3))

                if(haiopt.eq.1) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcicpx(1,1,1,4),qcicpy(1,1,1,4),        &
     &                          qcicpx_rst(1,1,1,4),qcicpy_rst(1,1,1,4))

                end if

              end if

            end if

! -----

! For the bin categories.

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qwtr(0,0,1,n),qwtr_rst(0,0,1,n))

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qwtrp(0,0,1,n),qwtrp_rst(0,0,1,n))

                if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            qwcpx(1,1,1,n),qwcpy(1,1,1,n),        &
     &                            qwcpx_rst(1,1,1,n),qwcpy_rst(1,1,1,n))

                end if

              end do

              do n=1,nnw

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          nwtr(0,0,1,n),nwtr_rst(0,0,1,n))

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          nwtrp(0,0,1,n),nwtrp_rst(0,0,1,n))

                if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            nwcpx(1,1,1,n),nwcpy(1,1,1,n),        &
     &                            nwcpx_rst(1,1,1,n),nwcpy_rst(1,1,1,n))

                end if

              end do

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        prwtr(0,0,1,1),prwtr_rst(0,0,1,1))

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qice(0,0,1,n),qice_rst(0,0,1,n))

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qicep(0,0,1,n),qicep_rst(0,0,1,n))

                if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            qicpx(1,1,1,n),qicpy(1,1,1,n),        &
     &                            qicpx_rst(1,1,1,n),qicpy_rst(1,1,1,n))

                end if

              end do

              do n=1,nni

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          nice(0,0,1,n),nice_rst(0,0,1,n))

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          nicep(0,0,1,n),nicep_rst(0,0,1,n))

                if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            nicpx(1,1,1,n),nicpy(1,1,1,n),        &
     &                            nicpx_rst(1,1,1,n),nicpy_rst(1,1,1,n))

                end if

              end do

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        price(0,0,1,1),price_rst(0,0,1,1))

            end if

          end if

! -----

!! -----

        end if

!!! -----

! For the aerosol variables.

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            call s_repsit3d('xx',istr,iend,jstr,jend,                   &
     &                      di,dj,ni,nj,nk,ni_rst,nj_rst,               &
     &                      qasl(0,0,1,n),qasl_rst(0,0,1,n))

            call s_repsit3d('xx',istr,iend,jstr,jend,                   &
     &                      di,dj,ni,nj,nk,ni_rst,nj_rst,               &
     &                      qaslp(0,0,1,n),qaslp_rst(0,0,1,n))

            if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

              call s_repsit2d('xx',istrb,iendb,jstrb,jendb,             &
     &                        dib,djb,ni,nj,nk,ni_rst,nj_rst,           &
     &                        qacpx(1,1,1,n),qacpy(1,1,1,n),            &
     &                        qacpx_rst(1,1,1,n),qacpy_rst(1,1,1,n))

            end if

          end do

        end if

! -----

! For the tracer variables.

        if(trkopt.ge.1) then

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,        &
     &                  ni_rst,nj_rst,qt,qt_rst)

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,        &
     &                  ni_rst,nj_rst,qtp,qtp_rst)

          if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

            call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,         &
     &                    ni,nj,nk,ni_rst,nj_rst,qtcpx,qtcpy,           &
     &                    qtcpx_rst,qtcpy_rst)

          end if

        end if

! -----

! For the turbulent kinetic energy variables.

        if(tubopt.ge.2) then

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,        &
     &                  ni_rst,nj_rst,tke,tke_rst)

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,        &
     &                  ni_rst,nj_rst,tkep,tkep_rst)

          if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

            call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,         &
     &                    ni,nj,nk,ni_rst,nj_rst,tkecpx,tkecpy,         &
     &                    tkecpx_rst,tkecpy_rst)

          end if

          if(dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-') then

            call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,      &
     &                    ni_rst,nj_rst,maxvl,maxvl_rst)

          end if

        end if

! -----

! For the diabatic variables.

        if(diaopt.eq.1) then

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,        &
     &                  ni_rst,nj_rst,pdia,pdia_rst)

        end if

! -----

! For the surface physics variables.

        if(sfcopt.ge.1) then

          call s_repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,1,       &
     &                    ni_rst,nj_rst,z0m,z0m_rst)

          call s_repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,1,       &
     &                    ni_rst,nj_rst,z0h,z0h_rst)

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nund,      &
     &                  ni_rst,nj_rst,tund,tund_rst)

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nund,      &
     &                  ni_rst,nj_rst,tundp,tundp_rst)

        end if

! -----

!!!! -----

!!!! Reposition the restructed variables in the case the Cubic Lagrange
!!!! advection scheme is performed.

      else

! For base state variables.

        call repsit3d('ox',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,ubr,ubr_rst)

        call repsit3d('xo',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,vbr,vbr_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,pbr,pbr_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,ptbr,ptbr_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,qvbr,qvbr_rst)

! -----

! For the velocity, the pressure and the potential temperature.

        call repsit3d('ox',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,up,up_rst)

        call repsit3d('xo',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,vp,vp_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,wp,wp_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,ppp,ppp_rst)

        call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,          &
     &                ni_rst,nj_rst,ptpp,ptpp_rst)

! -----

! For the boundary variables.

        if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

          call repsit2d('ox',istrb,iendb,jstrb,jendb,dib,djb,ni,nj,nk,  &
     &                  ni_rst,nj_rst,ucpx,ucpy,ucpx_rst,ucpy_rst)

          call repsit2d('xo',istrb,iendb,jstrb,jendb,dib,djb,ni,nj,nk,  &
     &                  ni_rst,nj_rst,vcpx,vcpy,vcpx_rst,vcpy_rst)

          call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,ni,nj,nk,  &
     &                  ni_rst,nj_rst,wcpx,wcpy,wcpx_rst,wcpy_rst)

          call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,ni,nj,nk,  &
     &                  ni_rst,nj_rst,pcpx,pcpy,pcpx_rst,pcpy_rst)

          call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,ni,nj,nk,  &
     &                  ni_rst,nj_rst,ptcpx,ptcpy,ptcpx_rst,ptcpy_rst)

        end if

! -----

!!! For the water substance.

        if(fmois(1:5).eq.'moist') then

! For the water vapor mixing ratio.

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,        &
     &                  ni_rst,nj_rst,qvp,qvp_rst)

          if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

            call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,         &
     &                    ni,nj,nk,ni_rst,nj_rst,qvcpx,qvcpy,           &
     &                    qvcpx_rst,qvcpy_rst)

          end if

! -----

!! For the water and ice hydrometeor.

! For the bulk categories.

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qwtrp(0,0,1,1),qwtrp_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qwtrp(0,0,1,2),qwtrp_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        prwtr(0,0,1,1),prwtr_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        prwtr(0,0,1,2),prwtr_rst(0,0,1,2))

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qwcpx(1,1,1,1),qwcpy(1,1,1,1),          &
     &                          qwcpx_rst(1,1,1,1),qwcpy_rst(1,1,1,1))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qwcpx(1,1,1,2),qwcpy(1,1,1,2),          &
     &                          qwcpx_rst(1,1,1,2),qwcpy_rst(1,1,1,2))

              end if

            end if

            if(abs(cphopt).eq.4) then

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nwtrp(0,0,1,1),nwtrp_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nwtrp(0,0,1,2),nwtrp_rst(0,0,1,2))

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nwcpx(1,1,1,1),nwcpy(1,1,1,1),          &
     &                          nwcpx_rst(1,1,1,1),nwcpy_rst(1,1,1,1))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nwcpx(1,1,1,2),nwcpy(1,1,1,2),          &
     &                          nwcpx_rst(1,1,1,2),nwcpy_rst(1,1,1,2))

              end if

            end if

            if(abs(cphopt).ge.2) then

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qicep(0,0,1,1),qicep_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qicep(0,0,1,2),qicep_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qicep(0,0,1,3),qicep_rst(0,0,1,3))

              if(haiopt.eq.1) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qicep(0,0,1,4),qicep_rst(0,0,1,4))

              end if

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        price(0,0,1,1),price_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        price(0,0,1,2),price_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        price(0,0,1,3),price_rst(0,0,1,3))

              if(haiopt.eq.1) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,2,ni_rst,nj_rst,            &
     &                          price(0,0,1,4),price_rst(0,0,1,4))

              end if

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qicpx(1,1,1,1),qicpy(1,1,1,1),          &
     &                          qicpx_rst(1,1,1,1),qicpy_rst(1,1,1,1))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qicpx(1,1,1,2),qicpy(1,1,1,2),          &
     &                          qicpx_rst(1,1,1,2),qicpy_rst(1,1,1,2))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qicpx(1,1,1,3),qicpy(1,1,1,3),          &
     &                          qicpx_rst(1,1,1,3),qicpy_rst(1,1,1,3))

                if(haiopt.eq.1) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            qicpx(1,1,1,4),qicpy(1,1,1,4),        &
     &                            qicpx_rst(1,1,1,4),qicpy_rst(1,1,1,4))

                end if

              end if

            end if

            if(abs(cphopt).eq.2) then

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nicep(0,0,1,1),nicep_rst(0,0,1,1))

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nicpx(1,1,1,1),nicpy(1,1,1,1),          &
     &                          nicpx_rst(1,1,1,1),nicpy_rst(1,1,1,1))

              end if

            else if(abs(cphopt).ge.3) then

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nicep(0,0,1,1),nicep_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nicep(0,0,1,2),nicep_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        nicep(0,0,1,3),nicep_rst(0,0,1,3))

              if(haiopt.eq.1) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          nicep(0,0,1,4),nicep_rst(0,0,1,4))

              end if

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nicpx(1,1,1,1),nicpy(1,1,1,1),          &
     &                          nicpx_rst(1,1,1,1),nicpy_rst(1,1,1,1))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nicpx(1,1,1,2),nicpy(1,1,1,2),          &
     &                          nicpx_rst(1,1,1,2),nicpy_rst(1,1,1,2))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          nicpx(1,1,1,3),nicpy(1,1,1,3),          &
     &                          nicpx_rst(1,1,1,3),nicpy_rst(1,1,1,3))

                if(haiopt.eq.1) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            nicpx(1,1,1,4),nicpy(1,1,1,4),        &
     &                            nicpx_rst(1,1,1,4),nicpy_rst(1,1,1,4))

                end if

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qcwtrp(0,0,1,1),qcwtrp_rst(0,0,1,1))

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qcwtrp(0,0,1,2),qcwtrp_rst(0,0,1,2))

              end if

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qcicep(0,0,1,1),qcicep_rst(0,0,1,1))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qcicep(0,0,1,2),qcicep_rst(0,0,1,2))

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,nk,ni_rst,nj_rst,             &
     &                        qcicep(0,0,1,3),qcicep_rst(0,0,1,3))

              if(haiopt.eq.1) then

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qcicep(0,0,1,4),qcicep_rst(0,0,1,4))

              end if

              if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                if(qcgopt.eq.2) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcwcpx(1,1,1,1),qcwcpy(1,1,1,1),        &
     &                          qcwcpx_rst(1,1,1,1),qcwcpy_rst(1,1,1,1))

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcwcpx(1,1,1,2),qcwcpy(1,1,1,2),        &
     &                          qcwcpx_rst(1,1,1,2),qcwcpy_rst(1,1,1,2))

                end if

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcicpx(1,1,1,1),qcicpy(1,1,1,1),        &
     &                          qcicpx_rst(1,1,1,1),qcicpy_rst(1,1,1,1))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcicpx(1,1,1,2),qcicpy(1,1,1,2),        &
     &                          qcicpx_rst(1,1,1,2),qcicpy_rst(1,1,1,2))

                call s_repsit2d('xx',istrb,iendb,jstrb,jendb,           &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcicpx(1,1,1,3),qcicpy(1,1,1,3),        &
     &                          qcicpx_rst(1,1,1,3),qcicpy_rst(1,1,1,3))

                if(haiopt.eq.1) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                          dib,djb,ni,nj,nk,ni_rst,nj_rst,         &
     &                          qcicpx(1,1,1,4),qcicpy(1,1,1,4),        &
     &                          qcicpx_rst(1,1,1,4),qcicpy_rst(1,1,1,4))

                end if

              end if

            end if

! -----

! For the bin categories.

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qwtrp(0,0,1,n),qwtrp_rst(0,0,1,n))

                if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            qwcpx(1,1,1,n),qwcpy(1,1,1,n),        &
     &                            qwcpx_rst(1,1,1,n),qwcpy_rst(1,1,1,n))

                end if

              end do

              do n=1,nnw

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          nwtrp(0,0,1,n),nwtrp_rst(0,0,1,n))

                if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            nwcpx(1,1,1,n),nwcpy(1,1,1,n),        &
     &                            nwcpx_rst(1,1,1,n),nwcpy_rst(1,1,1,n))

                end if

              end do

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        prwtr(0,0,1,1),prwtr_rst(0,0,1,1))

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          qicep(0,0,1,n),qicep_rst(0,0,1,n))

                if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            qicpx(1,1,1,n),qicpy(1,1,1,n),        &
     &                            qicpx_rst(1,1,1,n),qicpy_rst(1,1,1,n))

                end if

              end do

              do n=1,nni

                call s_repsit3d('xx',istr,iend,jstr,jend,               &
     &                          di,dj,ni,nj,nk,ni_rst,nj_rst,           &
     &                          nicep(0,0,1,n),nicep_rst(0,0,1,n))

                if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

                  call s_repsit2d('xx',istrb,iendb,jstrb,jendb,         &
     &                            dib,djb,ni,nj,nk,ni_rst,nj_rst,       &
     &                            nicpx(1,1,1,n),nicpy(1,1,1,n),        &
     &                            nicpx_rst(1,1,1,n),nicpy_rst(1,1,1,n))

                end if

              end do

              call s_repsit3d('xx',istr,iend,jstr,jend,                 &
     &                        di,dj,ni,nj,2,ni_rst,nj_rst,              &
     &                        price(0,0,1,1),price_rst(0,0,1,1))

            end if

          end if

! -----

!! -----

        end if

!!! -----

! For the aerosol variables.

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            call s_repsit3d('xx',istr,iend,jstr,jend,                   &
     &                      di,dj,ni,nj,nk,ni_rst,nj_rst,               &
     &                      qaslp(0,0,1,n),qaslp_rst(0,0,1,n))

            if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

              call s_repsit2d('xx',istrb,iendb,jstrb,jendb,             &
     &                        dib,djb,ni,nj,nk,ni_rst,nj_rst,           &
     &                        qacpx(1,1,1,n),qacpy(1,1,1,n),            &
     &                        qacpx_rst(1,1,1,n),qacpy_rst(1,1,1,n))

            end if

          end do

        end if

! -----

! For the tracer variables.

        if(trkopt.ge.1) then

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,        &
     &                  ni_rst,nj_rst,qtp,qtp_rst)

          if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

            call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,         &
     &                    ni,nj,nk,ni_rst,nj_rst,qtcpx,qtcpy,           &
     &                    qtcpx_rst,qtcpy_rst)

          end if

        end if

! -----

! For the turbulent kinetic energy variables.

        if(tubopt.ge.2) then

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,        &
     &                  ni_rst,nj_rst,tkep,tkep_rst)

          if(wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4) then

            call repsit2d('xx',istrb,iendb,jstrb,jendb,dib,djb,         &
     &                    ni,nj,nk,ni_rst,nj_rst,tkecpx,tkecpy,         &
     &                    tkecpx_rst,tkecpy_rst)

          end if

          if(dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-') then

            call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,      &
     &                    ni_rst,nj_rst,maxvl,maxvl_rst)

          end if

        end if

! -----

! For the diabatic variables.

        if(diaopt.eq.1) then

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nk,        &
     &                  ni_rst,nj_rst,pdia,pdia_rst)

        end if

! -----

! For the surface physics variables.

        if(sfcopt.ge.1) then

          call s_repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,1,       &
     &                    ni_rst,nj_rst,z0m,z0m_rst)

          call s_repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,1,       &
     &                    ni_rst,nj_rst,z0h,z0h_rst)

          call repsit3d('xx',istr,iend,jstr,jend,di,dj,ni,nj,nund,      &
     &                  ni_rst,nj_rst,tundp,tundp_rst)

        end if

! -----

      end if

!!!! -----

!!!!! -----

      end subroutine s_repdrv

!-----7--------------------------------------------------------------7--

      end module m_repdrv
