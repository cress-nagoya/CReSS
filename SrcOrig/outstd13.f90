!***********************************************************************
      module m_outstd13
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/01/15
!     Modification: 2002/07/15, 2003/05/19, 2004/05/31, 2004/06/10,
!                   2005/02/10, 2006/09/21, 2008/05/02, 2008/08/25,
!                   2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the messages to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd13, s_outstd13

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd13

        module procedure s_outstd13

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
      subroutine s_outstd13(sname,ncsn)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: sname
                       ! Called procedure name

      integer, intent(in) :: ncsn
                       ! Number of character of sname

!-----7--------------------------------------------------------------7--

! Read in the messages to standard i/o.

      write(6,*)

      write(6,'(a,a,a)') '  i/o: procedure, ',sname(1:ncsn),';'

      write(6,'(a)')                                                    &
     &        '    Read out user configuration from the standard input.'

! -----

      end subroutine s_outstd13

!-----7--------------------------------------------------------------7--

      end module m_outstd13
