!***********************************************************************
      module m_castname
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/04/06
!     Modification: 2000/01/17, 2003/05/19, 2008/05/02, 2008/08/25,
!                   2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     do nothing because this procedure is dummy.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: castname, s_castname

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface castname

        module procedure s_castname

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
      subroutine s_castname
!***********************************************************************

!-----7--------------------------------------------------------------7--

! Do nothing.

      end subroutine s_castname

!-----7--------------------------------------------------------------7--

      end module m_castname
