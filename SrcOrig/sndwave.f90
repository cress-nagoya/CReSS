!***********************************************************************
      module m_sndwave
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/04/06, 1999/07/05, 1999/08/18,
!                   1999/09/30, 1999/10/12, 1999/11/01, 2000/01/17,
!                   2000/07/05, 2002/04/02, 2003/01/04, 2003/04/30,
!                   2003/05/19, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the base state density x sound wave speed squared.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: sndwave, s_sndwave

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface sndwave

        module procedure s_sndwave

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
      subroutine s_sndwave(ni,nj,nk,pbr,rcsq)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

! Output variable

      real, intent(out) :: rcsq(0:ni+1,0:nj+1,1:nk)
                       ! rbr x sound wave speed squared

! Internal shared variable

      real cpdvcv      ! cp / cv

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Calculate cp / cv.

      cpdvcv=cp/cv

! ----

! Calculate the base state density x sound wave speed squared.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          rcsq(i,j,k)=cpdvcv*pbr(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_sndwave

!-----7--------------------------------------------------------------7--

      end module m_sndwave
