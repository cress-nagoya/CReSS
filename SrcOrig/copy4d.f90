!***********************************************************************
      module m_copy4d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/07/05, 1999/08/23, 1999/10/12,
!                   2000/01/17, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2006/01/10, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/08/18, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     copy 4 dimensional invar to outvar.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: copy4d, s_copy4d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface copy4d

        module procedure s_copy4d

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
      subroutine s_copy4d(imin,imax,jmin,jmax,kmin,kmax,nmin,nmax,      &
     &                    invar,outvar)
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

      integer, intent(in) :: nmin
                       ! Minimum array index in 4th direction

      integer, intent(in) :: nmax
                       ! Maximum array index in 4th direction

      real, intent(in) :: invar(imin:imax,jmin:jmax,kmin:kmax,nmin:nmax)
                       ! Coping variable

! Output variable

      real, intent(out) ::                                              &
     &                   outvar(imin:imax,jmin:jmax,kmin:kmax,nmin:nmax)
                       ! Copied variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in 4th direction

!-----7--------------------------------------------------------------7--

! Copy the invar to the outvar.

!$omp parallel default(shared) private(k,n)

      do n=nmin,nmax

        do k=kmin,kmax

!$omp do schedule(runtime) private(i,j)

          do j=jmin,jmax
          do i=imin,imax
            outvar(i,j,k,n)=invar(i,j,k,n)
          end do
          end do

!$omp end do

        end do

      end do

!$omp end parallel

! -----

      end subroutine s_copy4d

!-----7--------------------------------------------------------------7--

      end module m_copy4d
