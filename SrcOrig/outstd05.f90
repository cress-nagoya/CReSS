!***********************************************************************
      module m_outstd05
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/21
!     Modification: 1999/10/22, 2000/01/17, 2000/03/23, 2002/07/15,
!                   2003/05/19, 2004/05/31, 2004/06/10, 2004/08/20,
!                   2007/01/20, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in message to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comstd

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd05, s_outstd05

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd05

        module procedure s_outstd05

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
      subroutine s_outstd05(fmsg)
!***********************************************************************

! Input variable

      integer, intent(in) :: fmsg
                       ! Control flag of message type

!-----7--------------------------------------------------------------7--

! Read in messages to standard i/o.

      if(fstd(1:3).eq.'off') then

        if(fmsg.eq.2) then

          write(6,*)

        end if

      else if(fstd(1:3).eq.'act') then

        write(6,*)

        write(6,'(a)') '  #######################'

        if(fmsg.ge.1) then

          write(6,*)

        end if

        write(fstd(1:3),'(a3)') 'off'

      end if

! -----

      end subroutine s_outstd05

!-----7--------------------------------------------------------------7--

      end module m_outstd05
