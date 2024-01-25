!***********************************************************************
      module m_comcapt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/01/05
!     Modification: 2007/04/11, 2007/05/14, 2007/07/30, 2007/11/26,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/01/30,
!                   2009/02/27, 2009/08/20, 2011/06/01, 2011/08/18,
!                   2011/09/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the caption for dumped variables.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      public

! Exceptional access control

!     none

!-----7--------------------------------------------------------------7--

! Module variables

      character(len=60) capt(1:104)
                       ! Caption for dumped variable

      integer ncpt(1:104)
                       ! Number of character of capt

! Module procedure

!     none

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

!     none

!-----7--------------------------------------------------------------7--

      end module m_comcapt
