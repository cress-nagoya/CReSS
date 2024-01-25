!***********************************************************************
      module m_xy2ij
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/04/06, 1999/05/10, 1999/10/12,
!                   2000/01/17, 2002/04/02, 2002/07/03, 2003/04/30,
!                   2003/05/19, 2006/09/21, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the real indices from the 2 dimensional x and the y
!     coordinates at the model grid points.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getindx
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: xy2ij, s_xy2ij

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface xy2ij

        module procedure s_xy2ij

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
      subroutine s_xy2ij(fpdxiv,fpdyiv,xo,imin,imax,jmin,jmax,          &
     &                   x2d,y2d,ri,rj)
!***********************************************************************

! Input variables

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      real, intent(in) :: x2d(imin:imax,jmin:jmax)
                       ! 2 dimensional x coordinates

      real, intent(in) :: y2d(imin:imax,jmin:jmax)
                       ! 2 dimensional y coordinates

! Output variables

      real, intent(out) :: ri(imin:imax,jmin:jmax)
                       ! Real indices in x direction

      real, intent(out) :: rj(imin:imax,jmin:jmax)
                       ! Real indices in y direction

! Internal shared variables

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)

! -----

! Get the maximum and minimim indices of do loops.

      call getindx(xo,imin,imax,jmin,jmax,istr,iend,jstr,jend)

! -----

! Calculate the real indices from the 2 dimensional x and the y
! coordinates at the model grid points.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

      do j=jstr,jend
      do i=istr,iend
        ri(i,j)=x2d(i,j)*dxiv+2.e0
        rj(i,j)=y2d(i,j)*dyiv+2.e0
      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_xy2ij

!-----7--------------------------------------------------------------7--

      end module m_xy2ij
