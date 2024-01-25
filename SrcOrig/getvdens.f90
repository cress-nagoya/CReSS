!***********************************************************************
      module m_getvdens
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/04/01
!     Modification: 2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the inverse of base state density.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getvdens, s_getvdens

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getvdens

        module procedure s_getvdens

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
      subroutine s_getvdens(ni,nj,nk,rbr,rbv)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

! Output variable

      real, intent(out) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Calculate the inverse of base state density.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          rbv(i,j,k)=1.e0/rbr(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_getvdens

!-----7--------------------------------------------------------------7--

      end module m_getvdens
