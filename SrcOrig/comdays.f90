!***********************************************************************
      module m_comdays
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/11/05, 2004/10/12, 2006/03/06, 2006/08/08,
!                   2006/09/30, 2007/01/20, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the referenced number of days.

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

      integer ela(0:12)
                       ! Number of elapse of days from start of year

      integer elaitc(0:12)
                       ! Number of elapse of days
                       ! from start of intercalary year

      integer rem(1:12)
                       ! Number of remainder of days to end of year

      integer remitc(1:12)
                       ! Number of remainder of days
                       ! to end of intercalary year

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

      end module m_comdays
