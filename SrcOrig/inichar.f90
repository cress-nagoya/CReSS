!***********************************************************************
      module m_inichar
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/04/06, 2000/01/17, 2001/02/13, 2003/04/30,
!                   2003/05/19, 2006/09/21, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     initialize the optional character variable with space.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: inichar, s_inichar

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface inichar

        module procedure s_inichar

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
      subroutine s_inichar(opchar)
!***********************************************************************

! Output variable

      character(len=108), intent(out) :: opchar
                       ! Optional initialized character variable

!-----7--------------------------------------------------------------7--

! Fill in the optional character variable with space.

      write(opchar(1:54),'(a54)')                                       &
     &     '                                                      '

      write(opchar(55:108),'(a54)')                                     &
     &     '                                                      '

! -----

      end subroutine s_inichar

!-----7--------------------------------------------------------------7--

      end module m_inichar
