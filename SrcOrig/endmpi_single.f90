!***********************************************************************
      module m_endmpi
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2003/05/19, 2008/05/02, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the final message to the standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_outstd02

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: endmpi, s_endmpi

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface endmpi

        module procedure s_endmpi

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
      subroutine s_endmpi(stat)
!***********************************************************************

! Input variable

      integer, intent(in) :: stat
                       ! Runtime status

!-----7--------------------------------------------------------------7--

! Read in the final message to the standard i/o.

      call outstd02(stat)

! -----

      end subroutine s_endmpi

!-----7--------------------------------------------------------------7--

      end module m_endmpi
