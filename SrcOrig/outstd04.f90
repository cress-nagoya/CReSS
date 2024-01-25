!***********************************************************************
      module m_outstd04
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/14
!     Modification: 1999/06/21, 1999/10/22, 2000/01/05, 2000/01/17,
!                   2000/03/23, 2002/06/18, 2002/07/15, 2003/03/28,
!                   2003/05/19, 2004/05/31, 2004/06/10, 2004/08/20,
!                   2007/01/05, 2007/01/20, 2007/07/30, 2007/09/04,
!                   2008/05/02, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in message to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_comstd

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd04, s_outstd04

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd04

        module procedure s_outstd04

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_outstd04(fmsg,ctime)
!***********************************************************************

! Input variables

      integer, intent(in) :: fmsg
                       ! Control flag of message type

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

!-----7--------------------------------------------------------------7--

! Read in messages to standard i/o.

      if(fstd(1:3).eq.'off') then

        write(6,*)

        if(fmsg.eq.0) then

          write(6,'(a)') '  ##### Information #####'

        else if(fmsg.eq.1) then

          if(ctime.lt.10000_i8) then

            write(6,'(a,i1,a,i3.3,a)') '  ##### Information at time, ', &
     &                 ctime/1000_i8,'.',mod(ctime,1000_i8),' [s] #####'

          else if(ctime.lt.100000_i8) then

            write(6,'(a,i2,a,i3.3,a)') '  ##### Information at time, ', &
     &                 ctime/1000_i8,'.',mod(ctime,1000_i8),' [s] #####'

          else if(ctime.lt.1000000_i8) then

            write(6,'(a,i3,a,i3.3,a)') '  ##### Information at time, ', &
     &                 ctime/1000_i8,'.',mod(ctime,1000_i8),' [s] #####'

          else if(ctime.lt.10000000_i8) then

            write(6,'(a,i4,a,i3.3,a)') '  ##### Information at time, ', &
     &                 ctime/1000_i8,'.',mod(ctime,1000_i8),' [s] #####'

          else if(ctime.lt.100000000_i8) then

            write(6,'(a,i5,a,i3.3,a)') '  ##### Information at time, ', &
     &                 ctime/1000_i8,'.',mod(ctime,1000_i8),' [s] #####'

          else if(ctime.lt.1000000000_i8) then

            write(6,'(a,i6,a,i3.3,a)') '  ##### Information at time, ', &
     &                 ctime/1000_i8,'.',mod(ctime,1000_i8),' [s] #####'

          else if(ctime.lt.10000000000_i8) then

            write(6,'(a,i7,a,i3.3,a)') '  ##### Information at time, ', &
     &                 ctime/1000_i8,'.',mod(ctime,1000_i8),' [s] #####'

          else if(ctime.ge.10000000000_i8) then

            write(6,'(a,i8,a,i3.3,a)') '  ##### Information at time, ', &
     &                 ctime/1000_i8,'.',mod(ctime,1000_i8),' [s] #####'

          end if

        end if

        write(fstd(1:3),'(a3)') 'act'

      end if

! -----

      end subroutine s_outstd04

!-----7--------------------------------------------------------------7--

      end module m_outstd04
