!***********************************************************************
      module m_outstd08
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/01/15
!     Modification: 2001/04/15, 2001/05/29, 2001/12/10, 2002/04/02,
!                   2002/07/03, 2003/03/28, 2003/04/30, 2003/05/19,
!                   2003/09/01, 2003/12/12, 2004/03/05, 2004/06/10,
!                   2008/04/17, 2008/05/02, 2008/08/25, 2008/10/10,
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

      public :: outstd08, s_outstd08

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd08

        module procedure s_outstd08

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
      subroutine s_outstd08(fmsg)
!***********************************************************************

! Input variable

      integer, intent(in) :: fmsg
                       ! Control flag of message type

!-----7--------------------------------------------------------------7--

! Read in the messages to standard i/o.

      write(6,*)

      write(6,'(a)') '  messages: procedure, closedmp;'

      write(6,'(a)',advance='no')                                       &
     &            '    Because of wrong option of dmpvar,'

      if(fmsg.eq.1) then

        write(6,'(a)') ' no dumped variable.'

      else if(fmsg.eq.2) then

        write(6,'(a)') ' no 3d dumped variable.'

      else if(fmsg.eq.3) then

        write(6,'(a)') ' no dumped monitor variable.'

      end if

! -----

      end subroutine s_outstd08

!-----7--------------------------------------------------------------7--

      end module m_outstd08
