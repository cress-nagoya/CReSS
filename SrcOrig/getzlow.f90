!***********************************************************************
      module m_getzlow
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/04/01
!     Modification: 2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the z physical coordinates at lowest plane.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getzlow, s_getzlow

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getzlow

        module procedure s_getzlow

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
      subroutine s_getzlow(ni,nj,nk,zph,za)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

! Output variable

      real, intent(out) :: za(0:ni+1,0:nj+1)
                       ! z physical coordinates at lowest plane

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Calculate the z physical coordinates at lowest plane.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        za(i,j)=.5e0*(zph(i,j,3)-zph(i,j,2))
      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_getzlow

!-----7--------------------------------------------------------------7--

      end module m_getzlow
