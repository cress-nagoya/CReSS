!***********************************************************************
      module m_outstd02
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/06/14, 2000/01/17, 2002/06/18, 2003/05/19,
!                   2008/05/02, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the final message to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd02, s_outstd02

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd02

        module procedure s_outstd02

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
      subroutine s_outstd02(stat)
!***********************************************************************

! Input variable

      integer, intent(in) :: stat
                       ! Runtime status

!-----7--------------------------------------------------------------7--

! Read in the final message to standard i/o.

      write(6,*)

      if(stat.eq.0) then

        write(6,'(a)') 'This program stopped normally.'

      else

        write(6,'(a)') 'This program stopped abnormally.'

      end if

! -----

      end subroutine s_outstd02

!-----7--------------------------------------------------------------7--

      end module m_outstd02
