!***********************************************************************
      module m_culintg
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/04/03
!     Modification: 2006/05/12, 2006/06/21, 2006/07/21, 2006/11/06,
!                   2007/01/05, 2007/01/20, 2007/07/30, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2011/11/10

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the time steps integration
!     of the Cubic Lagrange advection.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_stepcul

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: culintg, s_culintg

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface culintg

        module procedure s_culintg

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
      subroutine s_culintg(fmois,nvstp,dtb,dtsep,                       &
     &                     ni,nj,nk,nqw,nnw,nqi,nni,nqa,                &
     &                     j31,j32,jcb8w,mf,mf8u,mf8v,uf,vf,wf,         &
     &                     ppf,ptpf,qvf,qwtrf,nwtrf,qicef,nicef,        &
     &                     qcwtrf,qcicef,qaslf,qtf,tkef,tmp1,tmp2,tmp3, &
     &                     tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: nvstp
                       ! Number of steps
                       ! of vertical Cubic Lagrange advection

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

      real, intent(in) :: dtsep
                       ! Time steps interval
                       ! of vertical Cubic Lagrange advection

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

! Input and output variables

      real, intent(inout) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(inout) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(inout) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

      real, intent(inout) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(inout) :: qwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at future

      real, intent(inout) :: nwtrf(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at future

      real, intent(inout) :: qicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at future

      real, intent(inout) :: nicef(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at future

      real, intent(inout) :: qcwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water hydrometeor
                       ! at future

      real, intent(inout) :: qcicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice hydrometeor
                       ! at future

      real, intent(inout) :: qaslf(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at future

      real, intent(inout) :: qtf(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at future

      real, intent(inout) :: tkef(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at future

! Internal shared variables

      integer ivstp    ! Index of time steps
                       ! of vertical Cubic Lagrange advection

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

      real, intent(inout) :: tmp8(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp9(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp10(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! The loop for the Cubic Lagrange advection.

      do ivstp=1,nvstp

        call stepcul(idcphopt,idhaiopt,idqcgopt,idaslopt,idtrkopt,      &
     &               idtubopt,fmois,nvstp,ivstp,dtb,dtsep,ni,nj,nk,     &
     &               nqw,nnw,nqi,nni,nqa,j31,j32,jcb8w,mf,mf8u,mf8v,    &
     &               uf,vf,wf,ppf,ptpf,qvf,qwtrf,nwtrf,qicef,nicef,     &
     &               qcwtrf,qcicef,qaslf,qtf,tkef,tmp1,tmp2,tmp3,       &
     &               tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10)

      end do

! -----

      end subroutine s_culintg

!-----7--------------------------------------------------------------7--

      end module m_culintg
