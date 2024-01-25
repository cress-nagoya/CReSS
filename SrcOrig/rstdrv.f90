!***********************************************************************
      module m_rstdrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/07/01
!     Modification: 2005/08/05, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/07/21, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/05/07, 2007/07/30, 2008/01/11, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2009/12/05, 2011/08/18, 2011/09/22, 2011/11/10,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     restruct the restart files.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_comkind
      use m_commpi
      use m_currpe
      use m_getiname
      use m_getrname
      use m_mvrst
      use m_outrst
      use m_outstd05
      use m_rdresrst
      use m_repdrv
      use m_rmres
      use m_rstindx
      use m_rststep

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rstdrv, s_rstdrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rstdrv

        module procedure s_rstdrv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rstdrv(fprmopt_rst,fpflitv_rst,ni,nj,nk,             &
     &               nqw,nnw,nqi,nni,nqa,nund,ni_rst,nj_rst,ubr,vbr,    &
     &               pbr,ptbr,qvbr,u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,      &
     &               qv,qvp,qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep,&
     &               qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp,       &
     &               tke,tkep,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy,  &
     &               ptcpx,ptcpy,qvcpx,qvcpy,qwcpx,qwcpy,nwcpx,nwcpy,   &
     &               qicpx,qicpy,nicpx,nicpy,qcwcpx,qcwcpy,             &
     &               qcicpx,qcicpy,qacpx,qacpy,qtcpx,qtcpy,             &
     &               tkecpx,tkecpy,maxvl,prwtr,price,pdia,z0m,z0h,      &
     &               tund,tundp,ubr_rst,vbr_rst,pbr_rst,ptbr_rst,       &
     &               qvbr_rst,u_rst,up_rst,v_rst,vp_rst,w_rst,wp_rst,   &
     &               pp_rst,ppp_rst,ptp_rst,ptpp_rst,qv_rst,qvp_rst,    &
     &               qwtr_rst,qwtrp_rst,nwtr_rst,nwtrp_rst,             &
     &               qice_rst,qicep_rst,nice_rst,nicep_rst,             &
     &               qcwtr_rst,qcwtrp_rst,qcice_rst,qcicep_rst,         &
     &               qasl_rst,qaslp_rst,qt_rst,qtp_rst,tke_rst,tkep_rst,&
     &               ucpx_rst,ucpy_rst,vcpx_rst,vcpy_rst,               &
     &               wcpx_rst,wcpy_rst,pcpx_rst,pcpy_rst,               &
     &               ptcpx_rst,ptcpy_rst,qvcpx_rst,qvcpy_rst,           &
     &               qwcpx_rst,qwcpy_rst,nwcpx_rst,nwcpy_rst,           &
     &               qicpx_rst,qicpy_rst,nicpx_rst,nicpy_rst,           &
     &               qcwcpx_rst,qcwcpy_rst,qcicpx_rst,qcicpy_rst,       &
     &               qacpx_rst,qacpy_rst,qtcpx_rst,qtcpy_rst,           &
     &               tkecpx_rst,tkecpy_rst,maxvl_rst,                   &
     &               prwtr_rst,price_rst,pdia_rst,z0m_rst,z0h_rst,      &
     &               tund_rst,tundp_rst)
!***********************************************************************

! Input variables

      integer, intent(in) :: fprmopt_rst
                       ! Formal parameter of unique index of rmopt_rst

      integer, intent(in) :: fpflitv_rst
                       ! Formal parameter of unique index of flitv_rst

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

! Internal shared variables

      character(len=5) fmois
                       ! Control flag of air moisture

      integer rmopt_rst
                       ! Option for original restart files removing

      integer ies      ! Start index in entire domain in x direction
      integer jes      ! Start index in entire domain in y direction

      integer ies_rst  ! Start index in entire domain in x direction
                       ! for restructed files

      integer jes_rst  ! Start index in entire domain in y direction
                       ! for restructed files

      integer imin     ! Minimum index in entire domain in x direction
      integer imax     ! Maximum index in entire domain in x direction
      integer jmin     ! Minimum index in entire domain in y direction
      integer jmax     ! Maximum index in entire domain in y direction

      integer imin_rst ! Minimum index in entire domain
                       ! in x direction for restructed files

      integer imax_rst ! Maximum index in entire domain
                       ! in x direction for restructed files

      integer jmin_rst ! Minimum index in entire domain
                       ! in y direction for restructed files

      integer jmax_rst ! Maximum index in entire domain
                       ! in y direction for restructed files

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction

      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      integer di       ! Differential index to istr
      integer dj       ! Differential index to jstr

      integer istrb    ! Minimum do loops index in x direction
                       ! of lateral boundary

      integer iendb    ! Maximum do loops index in x direction
                       ! of lateral boundary

      integer jstrb    ! Minimum do loops index in y direction
                       ! of lateral boundary

      integer jendb    ! Maximum do loops index in y direction
                       ! of lateral boundary

      integer dib      ! Differential index to istrb
      integer djb      ! Differential index to jstrb

      integer(kind=i8) it
                       ! Index of main do loop

      integer(kind=i8) nstp0
                       ! Start index of main do loop

      integer(kind=i8) nstp1
                       ! End index of main do loop

      integer(kind=i8) ctime
                       ! Model current forecast time

      integer(kind=i8) fl103
                       ! 1000 x int(flitv_rst + 0.1)

      real flitv_rst   ! Time interval of processed file

      real, intent(inout) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(inout) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(inout) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(inout) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(inout) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(inout) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(inout) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(inout) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(inout) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(inout) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(inout) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(inout) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at present

      real, intent(inout) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(inout) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(inout) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(inout) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(inout) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at present

      real, intent(inout) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(inout) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at present

      real, intent(inout) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(inout) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

      real, intent(inout) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(inout) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at present

      real, intent(inout) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(inout) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at present

      real, intent(inout) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(inout) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at present

      real, intent(inout) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(inout) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(inout) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(inout) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at present

      real, intent(inout) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

      real, intent(inout) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(inout) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(inout) :: ucpx(1:nj,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on west and east boundary

      real, intent(inout) :: ucpy(1:ni,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on south and north boundary

      real, intent(inout) :: vcpx(1:nj,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on west and east boundary

      real, intent(inout) :: vcpy(1:ni,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on south and north boundary

      real, intent(inout) :: wcpx(1:nj,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on west and east boundary

      real, intent(inout) :: wcpy(1:ni,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on south and north boundary

      real, intent(inout) :: pcpx(1:nj,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on west and east boundary

      real, intent(inout) :: pcpy(1:ni,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on south and north boundary

      real, intent(inout) :: ptcpx(1:nj,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on west and east boundary

      real, intent(inout) :: ptcpy(1:ni,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on south and north boundary

      real, intent(inout) :: qvcpx(1:nj,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on west and east boundary

      real, intent(inout) :: qvcpy(1:ni,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on south and north boundary

      real, intent(inout) :: qwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on west and east boundary

      real, intent(inout) :: qwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on south and north boundary

      real, intent(inout) :: nwcpx(1:nj,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on west and east boundary

      real, intent(inout) :: nwcpy(1:ni,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on south and north boundary

      real, intent(inout) :: qicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on west and east boundary

      real, intent(inout) :: qicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on south and north boundary

      real, intent(inout) :: nicpx(1:nj,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on west and east boundary

      real, intent(inout) :: nicpy(1:ni,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on south and north boundary

      real, intent(inout) :: qcwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution
                       ! for water on west and east boundary

      real, intent(inout) :: qcwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution
                       ! for water on south and north boundary

      real, intent(inout) :: qcicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution
                       ! for ice on west and east boundary

      real, intent(inout) :: qcicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution
                       ! for ice on south and north boundary

      real, intent(inout) :: qacpx(1:nj,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on west and east boundary

      real, intent(inout) :: qacpy(1:ni,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on south and north boundary

      real, intent(inout) :: qtcpx(1:nj,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on west and east boundary

      real, intent(inout) :: qtcpy(1:ni,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on south and north boundary

      real, intent(inout) :: tkecpx(1:nj,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on west and east boundary

      real, intent(inout) :: tkecpy(1:ni,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on south and north boundary

      real, intent(inout) :: maxvl(0:ni+1,0:nj+1,1:nk)
                       ! Maximum instantaneous wind velocity

      real, intent(inout) :: prwtr(0:ni+1,0:nj+1,1:2,1:nqw)
                       ! Precipitation and accumulation for water

      real, intent(inout) :: price(0:ni+1,0:nj+1,1:2,1:nqi)
                       ! Precipitation and accumulation for ice

      real, intent(inout) :: pdia(0:ni+1,0:nj+1,1:nk)
                       ! Diabatic value

      real, intent(inout) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(inout) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(inout) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at present

      real, intent(inout) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

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
                       ! nwtr in restructed domain

      real, intent(inout) :: nwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nnw)
                       ! nwtrp in restructed domain

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

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fprmopt_rst,rmopt_rst)
      call getrname(fpflitv_rst,flitv_rst)

! -----

! Read in the message to standard i/o.

      call outstd05(0)

! -----

! Calculate the number of steps of the main do loop.

      call rststep(idflitv_rst,idstime,idetime,nstp0,nstp1)

! -----

! Set the common used variable.

      fl103=1000_i8*int(flitv_rst+.1e0,i8)

! -----

!!!!!! The loop for the time integration.

      do it=nstp0,nstp1

! Calculate the current forecast time.

        ctime=fl103*(it-1_i8)

! -----

!!!!! The loop for the number of group domain.

        do mysrl=0,nsrl-1

! Set the parameters of parallelizing.

          call currpe('rstruct ',7,'ijgrp')

! -----

!!!! The loop for the number of restructed domain.

          do mysub_rst=0,nsub_rst-1

! Set the parameters of parallelizing.

            call currpe('rstruct ',7,'ijrst')

! -----

!!! The loop for the number of original restart files.

            do mysub=0,nsub-1

! Set the parameters of parallelizing.

              call currpe('rstruct ',7,'ijsub')

! -----

! Calculate the start indices in entire domain.

              ies=(ni-3)*(nisub*igrp+isub)
              jes=(nj-3)*(njsub*jgrp+jsub)

              ies_rst=(ni_rst-3)*(nisub_rst*igrp+isub_rst)
              jes_rst=(nj_rst-3)*(njsub_rst*jgrp+jsub_rst)

! -----

! Calculate the minimum and maximum indices in entire domain.

              imin=ies+2
              jmin=jes+2

              imax=ies+(ni-2)
              jmax=jes+(nj-2)

              imin_rst=ies_rst+2
              jmin_rst=jes_rst+2

              imax_rst=ies_rst+(ni_rst-2)
              jmax_rst=jes_rst+(nj_rst-2)

! -----

!! Open and read out the data from the restart file and reposition the
!! restructed variables form original restart variables.

              if(imin.le.imax_rst.and.imax.ge.imin_rst                  &
     &          .and.jmin.le.jmax_rst.and.jmax.ge.jmin_rst) then

! Calculate the minimum and maximum and differential do loops indices to
! reposition.

                call rstindx(ni,nj,ies,jes,ies_rst,jes_rst,             &
     &                       imin,imax,jmin,jmax,imin_rst,imax_rst,     &
     &                       jmin_rst,jmax_rst,istr,iend,jstr,jend,     &
     &                       di,dj,istrb,iendb,jstrb,jendb,dib,djb)

! -----

! Open and read out the data from the restart file.

                call rdresrst(idexprim,idcrsdir,iddmpvar,idncexp,       &
     &                     idnccrs,idwbc,idebc,idsbc,idnbc,idsfcopt,    &
     &                     idadvopt,idcphopt,idqcgopt,idaslopt,idtrkopt,&
     &                     idtubopt,iddiaopt,ctime,fmois,ni,nj,nk,      &
     &                     nqw,nnw,nqi,nni,nqa,nund,ubr,vbr,pbr,ptbr,   &
     &                     qvbr,u,up,v,vp,w,wp,pp,ppp,ptp,ptpp,qv,qvp,  &
     &                     qwtr,qwtrp,nwtr,nwtrp,qice,qicep,nice,nicep, &
     &                     qcwtr,qcwtrp,qcice,qcicep,qasl,qaslp,qt,qtp, &
     &                     tke,tkep,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,      &
     &                     pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,           &
     &                     qwcpx,qwcpy,nwcpx,nwcpy,qicpx,qicpy,         &
     &                     nicpx,nicpy,qcwcpx,qcwcpy,qcicpx,qcicpy,     &
     &                     qacpx,qacpy,qtcpx,qtcpy,tkecpx,tkecpy,maxvl, &
     &                     prwtr,price,pdia,z0m,z0h,tund,tundp)

! -----

! Reposition the restructed variables form original restart variables.

                call repdrv(iddmpvar,                                   &
     &                  idwbc,idebc,idsbc,idnbc,idsfcopt,idadvopt,      &
     &                  idcphopt,idhaiopt,idqcgopt,idaslopt,idtrkopt,   &
     &                  idtubopt,iddiaopt,fmois,istr,iend,jstr,jend,    &
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

! -----
              end if

!! -----

            end do

!!! -----

! Reset the parameters of parallelizing.

            call currpe('rstruct ',7,'ijrst')

! -----

! Open and read in the data to the restructed restart file.

            call outrst(idexprim,idcrsdir,iddmpvar,idncexp,idnccrs,     &
     &               idwbc,idebc,idsbc,idnbc,idsfcopt,idadvopt,idcphopt,&
     &               idqcgopt,idaslopt,idtrkopt,idtubopt,iddiaopt,      &
     &               fmois,ctime,nk,nqw,nnw,nqi,nni,nqa,nund,           &
     &               ni_rst,nj_rst,ubr_rst,vbr_rst,pbr_rst,ptbr_rst,    &
     &               qvbr_rst,u_rst,up_rst,v_rst,vp_rst,w_rst,wp_rst,   &
     &               pp_rst,ppp_rst,ptp_rst,ptpp_rst,qv_rst,qvp_rst,    &
     &               qwtr_rst,qwtrp_rst,nwtr_rst,nwtrp_rst,             &
     &               qice_rst,qicep_rst,nice_rst,nicep_rst,             &
     &               qcwtr_rst,qcwtrp_rst,qcice_rst,qcicep_rst,         &
     &               qasl_rst,qaslp_rst,qt_rst,qtp_rst,tke_rst,tkep_rst,&
     &               ucpx_rst,ucpy_rst,vcpx_rst,vcpy_rst,               &
     &               wcpx_rst,wcpy_rst,pcpx_rst,pcpy_rst,               &
     &               ptcpx_rst,ptcpy_rst,qvcpx_rst,qvcpy_rst,           &
     &               qwcpx_rst,qwcpy_rst,nwcpx_rst,nwcpy_rst,           &
     &               qicpx_rst,qicpy_rst,nicpx_rst,nicpy_rst,           &
     &               qcwcpx_rst,qcwcpy_rst,qcicpx_rst,qcicpy_rst,       &
     &               qacpx_rst,qacpy_rst,qtcpx_rst,qtcpy_rst,           &
     &               tkecpx_rst,tkecpy_rst,maxvl_rst,                   &
     &               prwtr_rst,price_rst,pdia_rst,z0m_rst,z0h_rst,      &
     &               tund_rst,tundp_rst)

! -----

          end do

!!!! -----

!!! Remove the original restart files and rename restructed restart
!!! files.

          if(rmopt_rst.eq.1) then

!! Remove the original restart files.

            do mysub=0,nsub-1

! Set the parameters of parallelizing.

              call currpe('rstruct ',7,'ijsub')

! -----

! Perform removing.

              call rmres(idexprim,idcrsdir,idncexp,idnccrs,ctime)

! -----

            end do

!! -----

!! Rename restructed restart files.

            do mysub_rst=0,nsub_rst-1

! Set the parameters of parallelizing.

              call currpe('rstruct ',7,'ijrst')

! -----

! Perform renaming.

              call mvrst(idexprim,idcrsdir,iddmpvar,idncexp,idnccrs,    &
     &               idwbc,idebc,idsbc,idnbc,idsfcopt,idadvopt,idcphopt,&
     &               idqcgopt,idaslopt,idtrkopt,idtubopt,iddiaopt,      &
     &               ctime,fmois,nk,nqw,nnw,nqi,nni,nqa,nund,           &
     &               ni_rst,nj_rst,ubr_rst,vbr_rst,pbr_rst,ptbr_rst,    &
     &               qvbr_rst,u_rst,up_rst,v_rst,vp_rst,w_rst,wp_rst,   &
     &               pp_rst,ppp_rst,ptp_rst,ptpp_rst,qv_rst,qvp_rst,    &
     &               qwtr_rst,qwtrp_rst,nwtr_rst,nwtrp_rst,             &
     &               qice_rst,qicep_rst,nice_rst,nicep_rst,             &
     &               qcwtr_rst,qcwtrp_rst,qcice_rst,qcicep_rst,         &
     &               qasl_rst,qaslp_rst,qt_rst,qtp_rst,tke_rst,tkep_rst,&
     &               ucpx_rst,ucpy_rst,vcpx_rst,vcpy_rst,               &
     &               wcpx_rst,wcpy_rst,pcpx_rst,pcpy_rst,               &
     &               ptcpx_rst,ptcpy_rst,qvcpx_rst,qvcpy_rst,           &
     &               qwcpx_rst,qwcpy_rst,nwcpx_rst,nwcpy_rst,           &
     &               qicpx_rst,qicpy_rst,nicpx_rst,nicpy_rst,           &
     &               qcwcpx_rst,qcwcpy_rst,qcicpx_rst,qcicpy_rst,       &
     &               qacpx_rst,qacpy_rst,qtcpx_rst,qtcpy_rst,           &
     &               tkecpx_rst,tkecpy_rst,maxvl_rst,                   &
     &               prwtr_rst,price_rst,pdia_rst,z0m_rst,z0h_rst,      &
     &               tund_rst,tundp_rst)

! -----

            end do

!! -----

          end if

!!! -----

        end do

!!!!! -----

! Read in the message to standard i/o.

        call outstd05(0)

! -----

      end do

!!!!!! -----

      end subroutine s_rstdrv

!-----7--------------------------------------------------------------7--

      end module m_rstdrv
