!***********************************************************************
      module m_copy1d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/07/05, 1999/08/23, 1999/10/12,
!                   2000/01/17, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2006/01/10, 2007/01/31, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     copy 1 dimensional invar to outvar.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: copy1d, s_copy1d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface copy1d

        module procedure s_copy1d

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
      subroutine s_copy1d(kmin,kmax,invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      real, intent(in) :: invar(kmin:kmax)
                       ! Coping variable

! Output variable

      real, intent(out) :: outvar(kmin:kmax)
                       ! Copied variable

! Internal private variable

      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Copy the invar to the outvar.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k)

      do k=kmin,kmax
        outvar(k)=invar(k)
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_copy1d

!-----7--------------------------------------------------------------7--

      end module m_copy1d
