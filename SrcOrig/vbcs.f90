!***********************************************************************
      module m_vbcs
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/20,
!                   1999/08/18, 1999/09/06, 1999/10/12, 1999/11/01,
!                   1999/12/06, 2000/01/17, 2000/03/17, 2001/09/13,
!                   2001/11/20, 2001/12/11, 2002/04/02, 2002/07/23,
!                   2002/08/15, 2003/04/30, 2003/05/19, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the vertical boundary condition for optional scalar variable.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vbcs, s_vbcs

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vbcs

        module procedure s_vbcs

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
      subroutine s_vbcs(ni,nj,nk,sf)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variable

      real, intent(inout) :: sf(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at future

! Internal shared variables

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

! -----

!! Set the bottom and top boundary conditions.

!$omp parallel default(shared)

! Set the bottom boundary conditions.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        sf(i,j,1)=sf(i,j,2)
      end do
      end do

!$omp end do

! -----

! Set the top boundary conditions.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        sf(i,j,nkm1)=sf(i,j,nkm2)
      end do
      end do

!$omp end do

! -----

!$omp end parallel

!! -----

      end subroutine s_vbcs

!-----7--------------------------------------------------------------7--

      end module m_vbcs
