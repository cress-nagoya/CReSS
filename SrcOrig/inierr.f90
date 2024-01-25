!***********************************************************************
      module m_inierr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/11/05
!     Modification: 2006/12/04, 2007/01/05, 2007/01/20, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     initialize the error list table.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comerr

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: inierr, s_inierr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface inierr

        module procedure s_inierr

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
      subroutine s_inierr
!***********************************************************************

! Internal shared variable

      integer ierr     ! Index of table to check namelist error

!-----7--------------------------------------------------------------7--

! Initialize the error list table.

      do ierr=1,nerr

        write(errlst(ierr)(1:14),'(a14)') '              '

      end do

! -----

      end subroutine s_inierr

!-----7--------------------------------------------------------------7--

      end module m_inierr
