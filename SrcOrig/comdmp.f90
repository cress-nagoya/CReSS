!***********************************************************************
      module m_comdmp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/11/19, 2000/01/17, 2000/04/18,
!                   2000/07/05, 2001/05/29, 2002/06/18, 2003/04/30,
!                   2003/05/19, 2006/09/21, 2006/12/04, 2007/01/20,
!                   2008/01/11, 2008/04/17, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/02/27, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the control variables of dumped file.

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

      character(len=108) fl3c
                       ! Name of dumped data checking file

      character(len=108) fl2c
                       ! Name of dumped data checking file
                       ! for monitor variables

      character(len=108) fl3d
                       ! Name of dumped file

      character(len=108) fl2d
                       ! Name of dumped file for monitor variables

      character(len=3) fdmp
                       ! Control flag of marked time for dumped file

      character(len=3) fmon
                       ! Control flag of marked time for dumped file
                       ! for monitor variables

      character(len=7) border
                       ! Identifier of endian

      integer nc3c     ! Number of character
                       ! of dumped data checking file name

      integer nc2c     ! Number of character
                       ! of dumped data checking file name
                       ! for monitor variables

      integer nc3d     ! Number of character of dumped file name

      integer nc2d     ! Number of character of dumped file name
                       ! for monitor variables

      integer io3c     ! Unit number of dumped data checking file

      integer io2c     ! Unit number of dumped data checking file
                       ! for monitor variables

      integer io3d     ! Unit number of dumped file

      integer io2d     ! Unit number of dumped file
                       ! for monitor variables

      integer rec3d    ! Current record number of dumped file

      integer rec2d    ! Current record number of dumped file
                       ! for monitor variables

      integer cnt3d    ! Dumped variables counter

      integer cnt2d    ! Dumped variables counter
                       ! for monitor variables

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

      end module m_comdmp
