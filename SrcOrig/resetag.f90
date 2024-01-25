!***********************************************************************
      module m_resetag
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/10/12
!     Modification: 2000/01/17, 2003/05/19, 2008/05/02, 2008/08/25,
!                   2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     reset the message tag.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: resetag, s_resetag

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface resetag

        module procedure s_resetag

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
      subroutine s_resetag
!***********************************************************************

!-----7--------------------------------------------------------------7--

! Reset the message tag.

      tag=0

! -----

      end subroutine s_resetag

!-----7--------------------------------------------------------------7--

      end module m_resetag
