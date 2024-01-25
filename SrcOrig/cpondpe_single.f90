!***********************************************************************
      module m_cpondpe
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/08/03
!     Modification: 2003/05/19, 2008/05/02, 2008/08/25, 2009/02/27

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

      public :: cpondpe, s_cpondpe

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface cpondpe

        module procedure s_cpondpe

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
      subroutine s_cpondpe
!***********************************************************************

!-----7--------------------------------------------------------------7--

! Do nothing.

      end subroutine s_cpondpe

!-----7--------------------------------------------------------------7--

      end module m_cpondpe
