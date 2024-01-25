!***********************************************************************
      module m_outstd15
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/02/10
!     Modification: 2005/08/05, 2006/09/21, 2007/01/05, 2007/01/20,
!                   2008/04/17, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27

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

      public :: outstd15, s_outstd15

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd15

        module procedure s_outstd15

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
      subroutine s_outstd15(vname,ncvn,fmsg1,fmsg2,dx)
!***********************************************************************

! Input variables

      character(len=6), intent(in) :: vname
                       ! Optional processed variable name

      integer, intent(in) :: ncvn
                       ! Number of character of vname

      integer, intent(in) :: fmsg1
                       ! Control flag of message type

      integer, intent(in) :: fmsg2
                       ! Control flag of message type

      real, intent(in) :: dx
                       ! Grid distance in x direction

!-----7--------------------------------------------------------------7--

! Read in messages to standard i/o.

      if(fmsg1.eq.1) then

        write(6,*)

      end if

      if(fmsg2.eq.1) then

        write(6,'(a,a,a,e13.6e2,a)')                                    &
     &     '    Reset the grid distance in x direction, ',vname(1:ncvn),&
     &     ' = ',dx,' [m].'

      else

        write(6,'(a,a,a,e13.6e2,a)')                                    &
     &     '    Reset the grid distance in x direction, ',vname(1:ncvn),&
     &     ' = ',dx,' [degree].'

      end if

! -----

      end subroutine s_outstd15

!-----7--------------------------------------------------------------7--

      end module m_outstd15
