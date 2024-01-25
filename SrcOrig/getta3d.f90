!***********************************************************************
      module m_getta3d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/09/30
!     Modification: 2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the air temperature.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getta3d, s_getta3d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getta3d

        module procedure s_getta3d

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
      subroutine s_getta3d(ni,nj,nk,ptbr,pi,ptp,t)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

! Output variable

      real, intent(out) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Calculate the air temperature.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          t(i,j,k)=(ptbr(i,j,k)+ptp(i,j,k))*pi(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_getta3d

!-----7--------------------------------------------------------------7--

      end module m_getta3d
