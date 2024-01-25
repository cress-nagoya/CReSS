!***********************************************************************
      module m_outstd03
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/14
!     Modification: 1999/06/21, 1999/10/22, 2000/01/05, 2000/01/17,
!                   2000/03/23, 2000/07/05, 2001/05/29, 2002/06/18,
!                   2002/07/15, 2003/04/30, 2003/05/19, 2004/05/31,
!                   2004/06/10, 2005/02/10, 2006/09/21, 2007/07/30,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2013/03/27

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

      public :: outstd03, s_outstd03

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd03

        module procedure s_outstd03

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
      subroutine s_outstd03(sname,ncsn,fname,ncfn,ionum,fmsg1,fmsg2,    &
     &                      ctime)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: sname
                       ! Called procedure name

      character(len=108), intent(in) :: fname
                       ! Opned file name

      integer, intent(in) :: ncsn
                       ! Number of character of sname

      integer, intent(in) :: ncfn
                       ! Number of character of fname

      integer, intent(in) :: ionum
                       ! File unit number

      integer, intent(in) :: fmsg1
                       ! Control flag of message type

      integer, intent(in) :: fmsg2
                       ! Control flag of message type

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

!-----7--------------------------------------------------------------7--

! Read in the messages to standard i/o.

      call outstd04(fmsg2,ctime)

      write(6,*)

      write(6,'(a,a,a)') '  i/o: procedure, ',sname(1:ncsn),';'

      if(fmsg1.eq.1) then

        if(ncfn.le.32) then

         write(6,'(a,a,a,i5.5,a)') '    Open the file(s), "',           &
     &               fname(1:ncfn),'" with unit number ',ionum,'.'

        else

         write(6,'(a,a,a)') '    Open the file(s), "',fname(1:ncfn),'"'

         write(6,'(a,i5.5,a)') '    with unit number ',ionum,'.'

        end if

      else if(fmsg1.eq.2) then

        write(6,'(a,i5.5,a)')                                           &
     &       '    Closed the file(s) of unit number ',ionum,'.'

      else if(fmsg1.eq.3) then

        write(6,'(a,i5.5,a)')                                           &
     &       '    Read out the data from the file(s) of unit number ',  &
     &       ionum,'.'

      else if(fmsg1.eq.4) then

        write(6,'(a,i5.5,a)')                                           &
     &       '    Read in the data to the file(s) of unit number ',     &
     &       ionum,'.'

      end if

! -----

      end subroutine s_outstd03

!-----7--------------------------------------------------------------7--

      end module m_outstd03
