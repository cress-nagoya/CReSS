!***********************************************************************
      module m_getcname
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/03/25
!     Modification: 1999/04/06, 2000/01/17, 2002/04/02, 2003/04/30,
!                   2003/05/19, 2004/01/09, 2006/09/21, 2007/01/20,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the character namelist variable.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getcname, s_getcname

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getcname

        module procedure s_getcname

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
      subroutine s_getcname(idc,cval)
!***********************************************************************

! Input variable

      integer, intent(in) :: idc
                       ! Optional index in character namelist table

! Output variable

      character(len=108), intent(out) :: cval
                       ! Character namelist variable

!-----7--------------------------------------------------------------7--

! Get the character namelist variable.

      cval(1:108)=cname(idc)(1:108)

! -----

      end subroutine s_getcname

!-----7--------------------------------------------------------------7--

      end module m_getcname
