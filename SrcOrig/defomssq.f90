!***********************************************************************
      module m_defomssq
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/08/07
!     Modification: 2001/11/20, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2003/11/05, 2004/02/01, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the deformation squared.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: defomssq, s_defomssq

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface defomssq

        module procedure s_defomssq

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
      subroutine s_defomssq(ni,nj,nk,s11,s22,s33,s12,s31,s32,ssq)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: s11(0:ni+1,0:nj+1,1:nk)
                       ! x-x components of deformation tensor

      real, intent(in) :: s22(0:ni+1,0:nj+1,1:nk)
                       ! y-y components of deformation tensor

      real, intent(in) :: s33(0:ni+1,0:nj+1,1:nk)
                       ! z-z components of deformation tensor

      real, intent(in) :: s12(0:ni+1,0:nj+1,1:nk)
                       ! x-y components of deformation tensor

      real, intent(in) :: s31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of deformation tensor

      real, intent(in) :: s32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of deformation tensor

! Output variable

      real, intent(out) :: ssq(0:ni+1,0:nj+1,1:nk)
                       ! Magnitude of deformation squared

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real s128s       ! s12 at scalar points
      real s328s       ! s32 at scalar points
      real s318s       ! s31 at scalar points

!-----7--------------------------------------------------------------7--

! Calculate the magnitude of the deformation squared.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j,s128s,s318s,s328s)

        do j=1,nj-1
        do i=1,ni-1
          s128s=(s12(i,j,k)+s12(i+1,j+1,k))+(s12(i+1,j,k)+s12(i,j+1,k))
          s318s=(s31(i,j,k)+s31(i+1,j,k+1))+(s31(i,j,k+1)+s31(i+1,j,k))
          s328s=(s32(i,j,k)+s32(i,j+1,k+1))+(s32(i,j+1,k)+s32(i,j,k+1))

          ssq(i,j,k)=.5e0*(s33(i,j,k)*s33(i,j,k)                        &
     &      +(s11(i,j,k)*s11(i,j,k)+s22(i,j,k)*s22(i,j,k)))             &
     &      +.0625e0*(s128s*s128s+(s318s*s318s+s328s*s328s))

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_defomssq

!-----7--------------------------------------------------------------7--

      end module m_defomssq
