!***********************************************************************
      module m_bcphi
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/06/01
!     Modification: 2001/08/07, 2001/12/11, 2002/04/02, 2002/10/15,
!                   2003/04/30, 2003/05/19, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions for the parabolic pertial differential
!     equation.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bcphi, s_bcphi

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcphi

        module procedure s_bcphi

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
      subroutine s_bcphi(ni,nj,nk,phi)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Input and output variable

      real, intent(inout) :: phi(0:ni+1,0:nj+1,1:nk)
                       ! Optional solved variable

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

! Set the bottom and top boundary conditions.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        phi(i,j,1)=phi(i,j,2)
        phi(i,j,nkm1)=phi(i,j,nkm2)
      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_bcphi

!-----7--------------------------------------------------------------7--

      end module m_bcphi
