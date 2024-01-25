!***********************************************************************
      program check
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/05/31
!     Modification: 2005/02/10, 2007/01/20, 2007/07/30, 2008/04/17,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! This is the main program for the pre processor ckeck.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkkind
      use m_chkname
      use m_ini0name
      use m_inierr
      use m_inimpi
      use m_outstd02
      use m_outstd05

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Initialize the parameters of parallelizing.

      call inimpi('check   ',5)

! -----

! Check the kind of the Fortran variables.

      call chkkind

! -----

! Initialize the namelist variables.

      call ini0name

! -----

! Initialize the error list table.

      call inierr

! -----

! Read out and check the namelist variables.

      call chkname

! -----

! Read in the final message to standard i/o.

      call outstd05(0)
      call outstd02(0)

! -----

!-----7--------------------------------------------------------------7--

! Internal procedure

!     none

!-----7--------------------------------------------------------------7--

      end program check
