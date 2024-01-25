!***********************************************************************
      module m_getz11
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/03/25
!     Modification: 1999/04/06, 1999/05/10, 1999/11/01, 2000/01/17,
!                   2002/04/02, 2003/04/30, 2003/05/19, 2005/08/05,
!                   2007/01/31, 2007/06/27, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the zeta coordinates.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getz11, s_getz11

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getz11

        module procedure s_getz11

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getz11(fpdz,fpzsfc,nk,z)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: fpzsfc
                       ! Formal parameter of unique index of zsfc

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Output variable

      real, intent(out) :: z(1:nk)
                       ! zeta coordinates

! Internal shared variables

      real dz          ! Grid distance in z direction

      real zsfc        ! Sea surface terrain height

      real zsfc11      ! zsfc + 11.0

! Internal private variable

      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpdz,dz)
      call getrname(fpzsfc,zsfc)

! -----

! Set the common used variable.

      zsfc11=zsfc+11.e0

! -----

! Calculate the zeta coordinates.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k)

      do k=1,nk
        z(k)=zsfc11+real(k-2)*dz
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_getz11

!-----7--------------------------------------------------------------7--

      end module m_getz11
