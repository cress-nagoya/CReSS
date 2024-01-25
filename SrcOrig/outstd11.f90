!***********************************************************************
      module m_outstd11
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/08/07
!     Modification: 2002/06/18, 2002/07/15, 2002/10/15, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2004/05/31, 2004/06/10,
!                   2006/09/21, 2007/01/05, 2007/07/30, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the message to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_outstd04

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd11, s_outstd11

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd11

        module procedure s_outstd11

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
      subroutine s_outstd11(sname,ncsn,vcap,ncvc,itcnt,fmsg1,fmsg2,     &
     &                      ctime)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: sname
                       ! Called procedure name

      character(len=26), intent(in) :: vcap
                       ! Caption for processed variable

      integer, intent(in) :: ncsn
                       ! Number of character of sname

      integer, intent(in) :: ncvc
                       ! Number of character of vcap

      integer, intent(in) :: itcnt
                       ! Iteration count

      integer, intent(in) :: fmsg1
                       ! Control flag of message type

      integer, intent(in) :: fmsg2
                       ! Control flag of message type

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

!-----7--------------------------------------------------------------7--

! Read in the messages to standart i/o.

      call outstd04(fmsg2,ctime)

      write(6,*)

      write(6,'(a,a,a)') '  iteration: procedure, ',sname(1:ncsn),';'

      if(fmsg1.eq.1) then

        write(6,'(a,a,a)') '    Start ',vcap(1:ncvc),' fitting.'

      else if(fmsg1.eq.2) then

        if(itcnt.lt.10) then

          write(6,'(a,a,a,i1,a)') '    End ',vcap(1:ncvc),              &
     &                            ' fitting, iterated ',itcnt,' times.'

        else if(itcnt.lt.100) then

          write(6,'(a,a,a,i2,a)') '    End ',vcap(1:ncvc),              &
     &                            ' fitting, iterated ',itcnt,' times.'

        else if(itcnt.lt.1000) then

          write(6,'(a,a,a,i3,a)') '    End ',vcap(1:ncvc),              &
     &                            ' fitting, iterated ',itcnt,' times.'

        else if(itcnt.lt.10000) then

          write(6,'(a,a,a,i4,a)') '    End ',vcap(1:ncvc),              &
     &                            ' fitting, iterated ',itcnt,' times.'

        else if(itcnt.lt.100000) then

          write(6,'(a,a,a,i5,a)') '    End ',vcap(1:ncvc),              &
     &                            ' fitting, iterated ',itcnt,' times.'

        else if(itcnt.lt.1000000) then

          write(6,'(a,a,a,i6,a)') '    End ',vcap(1:ncvc),              &
     &                            ' fitting, iterated ',itcnt,' times.'

        end if

      end if

! -----

      end subroutine s_outstd11

!-----7--------------------------------------------------------------7--

      end module m_outstd11
