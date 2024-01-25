!***********************************************************************
      module m_comfile
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2008/10/10
!     Modification: 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     declare the parameters of constant directory and file name.

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

      character(len=2), parameter :: curdir='./'
                       ! Current directory name

      character(len=11), parameter :: klfl='kill.solver'
                       ! File name to abort solver

      character(len=16), parameter :: bytefl='endian-check.bin'
                       ! File name to check endian

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

      end module m_comfile
