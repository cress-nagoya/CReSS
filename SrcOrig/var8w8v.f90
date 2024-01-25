!***********************************************************************
      module m_var8w8v
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/05/07
!     Modification: 2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     optional variable at w points be averaged to v points.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: var8w8v, s_var8w8v

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface var8w8v

        module procedure s_var8w8v

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
      subroutine s_var8w8v(ni,nj,nk,var8w,var8v)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: var8w(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at w points

! Output variable

      real, intent(out) :: var8v(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at v points

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Be averaged to v points.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj
        do i=0,ni
          var8v(i,j,k)=.25e0*((var8w(i,j-1,k)+var8w(i,j-1,k+1))         &
     &      +(var8w(i,j,k)+var8w(i,j,k+1)))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_var8w8v

!-----7--------------------------------------------------------------7--

      end module m_var8w8v
