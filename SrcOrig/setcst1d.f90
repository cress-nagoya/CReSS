!***********************************************************************
      module m_setcst1d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/07/05, 1999/08/23, 1999/09/30,
!                   1999/10/12, 2000/01/17, 2002/04/02, 2003/04/30,
!                   2003/05/19, 2006/01/10, 2007/01/31, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2010/12/13,
!                   2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     fill in the 1 dimensional array with the constant value.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setcst1d, s_setcst1d, s_setcst1d_r8

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setcst1d

        module procedure s_setcst1d, s_setcst1d_r8

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
      subroutine s_setcst1d(kmin,kmax,invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      real, intent(in) :: invar
                       ! Constant variable

! Output variable

      real, intent(out) :: outvar(kmin:kmax)
                       ! Copied variable

! Internal private variable

      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Fill in the array with the constant value.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k)

      do k=kmin,kmax
        outvar(k)=invar
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_setcst1d

!***********************************************************************
      subroutine s_setcst1d_r8(kmin,kmax,invar,outvar)
!***********************************************************************

! Input variables

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      real, intent(in) :: invar
                       ! Constant variable

! Output variable

      real(kind=r8), intent(out) :: outvar(kmin:kmax)
                       ! Copied variable

! Internal private variable

      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Fill in the array with the constant value.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k)

      do k=kmin,kmax
        outvar(k)=invar
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_setcst1d_r8

!-----7--------------------------------------------------------------7--

      end module m_setcst1d
