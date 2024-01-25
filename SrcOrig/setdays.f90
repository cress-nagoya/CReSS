!***********************************************************************
      module m_setdays
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/11/05
!     Modification: 2004/10/12, 2004/12/17, 2006/03/06, 2006/08/08,
!                   2006/09/30, 2007/01/20, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the referenced number of days.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comdays

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setdays, s_setdays

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setdays

        module procedure s_setdays

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
      subroutine s_setdays
!***********************************************************************

!-----7--------------------------------------------------------------7--

! Set the table for the number of elapse of days form the start of year.

      ela(0)=0
      ela(1)=31
      ela(2)=59
      ela(3)=90
      ela(4)=120
      ela(5)=151
      ela(6)=181
      ela(7)=212
      ela(8)=243
      ela(9)=273
      ela(10)=304
      ela(11)=334
      ela(12)=365

! -----

! Set the table for the number of elapse of days form the start of
! intercalary year.

      elaitc(0)=0
      elaitc(1)=31
      elaitc(2)=60
      elaitc(3)=91
      elaitc(4)=121
      elaitc(5)=152
      elaitc(6)=182
      elaitc(7)=213
      elaitc(8)=244
      elaitc(9)=274
      elaitc(10)=305
      elaitc(11)=335
      elaitc(12)=366

! -----

! Set the table for the number of remainder of days to the end of year.

      rem(1)=365
      rem(2)=334
      rem(3)=306
      rem(4)=275
      rem(5)=245
      rem(6)=214
      rem(7)=184
      rem(8)=153
      rem(9)=122
      rem(10)=92
      rem(11)=61
      rem(12)=31

! -----

! Set the table for the number of remainder of days to the end of
! intercalary year.

      remitc(1)=366
      remitc(2)=335
      remitc(3)=306
      remitc(4)=275
      remitc(5)=245
      remitc(6)=214
      remitc(7)=184
      remitc(8)=153
      remitc(9)=122
      remitc(10)=92
      remitc(11)=61
      remitc(12)=31

! -----

      end subroutine s_setdays

!-----7--------------------------------------------------------------7--

      end module m_setdays
