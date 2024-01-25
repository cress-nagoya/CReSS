!***********************************************************************
      module m_outstd16
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/05/31
!     Modification: 2004/06/10, 2005/02/10, 2006/09/21, 2006/12/04,
!                   2007/01/20, 2008/04/17, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in message to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd16, s_outstd16

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd16

        module procedure s_outstd16

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
      subroutine s_outstd16(pname,ncpn,stat)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: pname
                       ! Running program name

      integer, intent(in) :: ncpn
                       ! Number of character of pname

      integer, intent(in) :: stat
                       ! Runtime status

!-----7--------------------------------------------------------------7--

! Read in messages to standard i/o.

      write(6,*)

      write(6,'(a)') '  messages: procedure, chkname;'

      if(stat.eq.0) then

       write(6,'(a,a,a,a)') '    There is no error in the',             &
     &         ' configuration file for the program, ',pname(1:ncpn),'.'

      else if(stat.eq.1) then

       write(6,'(a,i1,a,a,a,a)') '    There is ',stat,' error in the',  &
     &         ' configuration file for the program, ',pname(1:ncpn),'.'

      else if(stat.ge.2.and.stat.lt.10) then

       write(6,'(a,i1,a,a,a,a)') '    There are ',stat,' errors in the',&
     &         ' configuration file for the program, ',pname(1:ncpn),'.'

      else if(stat.ge.10.and.stat.lt.100) then

       write(6,'(a,i2,a,a,a,a)') '    There are ',stat,' errors in the',&
     &         ' configuration file for the program, ',pname(1:ncpn),'.'

      else

       write(6,'(a,i3,a,a,a,a)') '    There are ',stat,' errors in the',&
     &         ' configuration file for the program, ',pname(1:ncpn),'.'

      end if

! -----

      end subroutine s_outstd16

!-----7--------------------------------------------------------------7--

      end module m_outstd16
