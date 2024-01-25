!***********************************************************************
      module m_var8v8s
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/08/18
!     Modification: 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     optional variable at v points be averaged to scalar points.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: var8v8s, s_var8v8s

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface var8v8s

        module procedure s_var8v8s

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
      subroutine s_var8v8s(ni,nj,nk,var8v,var8s)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: var8v(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at v points

! Output variable

      real, intent(out) :: var8s(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at scalar points

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Be averaged to scalar points.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          var8s(i,j,k)=.5e0*(var8v(i,j,k)+var8v(i,j+1,k))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_var8v8s

!-----7--------------------------------------------------------------7--

      end module m_var8v8s
