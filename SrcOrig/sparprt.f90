!***********************************************************************
      module m_sparprt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/07/05, 1999/09/30, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2001/01/15, 2001/06/29, 2002/04/02,
!                   2003/04/30, 2003/05/19, 2004/01/09, 2006/02/03,
!                   2007/10/19, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2009/02/27, 2011/12/17, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     separate the pressure and the potential temperature to the base
!     state and the perturbation value.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_vint13

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: sparprt, s_sparprt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface sparprt

        module procedure s_sparprt

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
      subroutine s_sparprt(nid,njd,nkd,zdat,ppdat,ptpdat,pbdat,ptbdat,  &
     &                     nk,z1d,p1d,pt1d)
!***********************************************************************

! Input variables

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(in) :: z1d(0:4*nk-3)
                       ! Horizontally averaged z coordinates

      real, intent(in) :: p1d(0:4*nk-3)
                       ! Horizontally averaged pressure

      real, intent(in) :: pt1d(0:4*nk-3)
                       ! Horizontally averaged potential temperature

! Input and output variables

      real, intent(inout) :: ppdat(1:nid,1:njd,1:nkd)
                       ! Pressure perturbation in data

      real, intent(inout) :: ptpdat(1:nid,1:njd,1:nkd)
                       ! Potential temperature perturbation in data

! Internal shared variables

      integer nk4m3    ! 4 x nk - 3

      real, intent(inout) :: pbdat(1:nid,1:njd,1:nkd)
                       ! Base state pressure in data

      real, intent(inout) :: ptbdat(1:nid,1:njd,1:nkd)
                       ! Base state potential temperature in data

! Internal private variables

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction
      integer kd       ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Set the common used variable.

      nk4m3=4*nk-3

! -----

! Get the base state pressure and potential temperature in data grid.

      call vint13('oo',1,nid,1,njd,1,nkd,zdat,pbdat,nk4m3,z1d,p1d)

      call vint13('oo',1,nid,1,njd,1,nkd,zdat,ptbdat,nk4m3,z1d,pt1d)

! -----

! Separate the pressure and the potential temperature to the base state
! and the perturbation value.

!$omp parallel default(shared) private(kd)

      do kd=1,nkd

!$omp do schedule(runtime) private(id,jd)

        do jd=1,njd
        do id=1,nid
          ppdat(id,jd,kd)=ppdat(id,jd,kd)-pbdat(id,jd,kd)
          ptpdat(id,jd,kd)=ptpdat(id,jd,kd)-ptbdat(id,jd,kd)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_sparprt

!-----7--------------------------------------------------------------7--

      end module m_sparprt
