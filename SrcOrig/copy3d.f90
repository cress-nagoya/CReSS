!***********************************************************************
      module m_copy3d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/07/05, 1999/08/23, 1999/10/12,
!                   2000/01/17, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2006/01/10, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     copy 3 dimensional invar to outvar.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: copy3d, s_copy3d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface copy3d

        module procedure s_copy3d

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
      subroutine s_copy3d(imin,imax,jmin,jmax,kmin,kmax,invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      real, intent(in) :: invar(imin:imax,jmin:jmax,kmin:kmax)
                       ! Coping variable

! Output variable

      real, intent(out) :: outvar(imin:imax,jmin:jmax,kmin:kmax)
                       ! Copied variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Copy the invar to the outvar.

!$omp parallel default(shared) private(k)

      do k=kmin,kmax

!$omp do schedule(runtime) private(i,j)

        do j=jmin,jmax
        do i=imin,imax
          outvar(i,j,k)=invar(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_copy3d

!-----7--------------------------------------------------------------7--

      end module m_copy3d
