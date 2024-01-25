!***********************************************************************
      module m_vsps0
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/07/05, 1999/08/03, 1999/09/30, 1999/10/12,
!                   1999/11/01, 2000/01/17, 2001/04/15, 2002/04/02,
!                   2002/08/15, 2003/04/30, 2003/05/19, 2003/06/27,
!                   2003/12/12, 2004/04/15, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the vertical sponge damping for optional scalar variable
!     to initial.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vsps0, s_vsps0

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vsps0

        module procedure s_vsps0

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
      subroutine s_vsps0(ksp0,ni,nj,nk,rst,sp,rbct,sfrc)
!***********************************************************************

! Input variables

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: sp(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at past

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

! Input and output variable

      real, intent(inout) :: sfrc(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar forcing term

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Calculate the vertical sponge damping for optional scalar variable to
! initial.

!$omp parallel default(shared) private(k)

      do k=ksp0(2)-1,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          sfrc(i,j,k)=sfrc(i,j,k)                                       &
     &      -.5e0*(rbct(i,j,k,2)+rbct(i,j,k+1,2))*rst(i,j,k)*sp(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_vsps0

!-----7--------------------------------------------------------------7--

      end module m_vsps0
