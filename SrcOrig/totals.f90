!***********************************************************************
      module m_totals
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/08/18
!     Modification: 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the total value of optional scalar variable from base state
!     and perturbation value.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: totals, s_totals

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface totals

        module procedure s_totals

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
      subroutine s_totals(ni,nj,nk,sbr,sp,s)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: sbr(0:ni+1,0:nj+1,1:nk)
                       ! Optional base state scalar variable

      real, intent(in) :: sp(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar perturbation variable

! Output variable

      real, intent(out) :: s(0:ni+1,0:nj+1,1:nk)
                       ! Total value of optional scalar variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the total value of optional scalar variable from base state and
! perturbation value.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          s(i,j,k)=sbr(i,j,k)+sp(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_totals

!-----7--------------------------------------------------------------7--

      end module m_totals
