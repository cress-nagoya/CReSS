!***********************************************************************
      module m_chkstd
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/06/18
!     Modification: 2003/05/19, 2004/05/31, 2007/01/20, 2008/05/02,
!                   2008/08/25, 2009/02/27

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

      public :: chkstd, s_chkstd

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkstd

        module procedure s_chkstd

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
      subroutine s_chkstd(broot)
!***********************************************************************

! Input and output variable

      integer, intent(inout) :: broot
                       ! Broadcasting root

!-----7--------------------------------------------------------------7--

! Return same value.

      broot=1*broot

! -----

      end subroutine s_chkstd

!-----7--------------------------------------------------------------7--

      end module m_chkstd
