!***********************************************************************
      module m_getiname
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/03/25
!     Modification: 2000/01/17, 2002/04/02, 2003/05/19, 2004/01/09,
!                   2007/01/20, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the integer namelist variable.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getiname, s_getiname

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getiname

        module procedure s_getiname

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
      subroutine s_getiname(idi,ival)
!***********************************************************************

! Input variable

      integer, intent(in) :: idi
                       ! Optional index in integer namelist table

! Output variable

      integer, intent(out) :: ival
                       ! Integer namelist variable

!-----7--------------------------------------------------------------7--

! Get the integer namelist variable.

      ival=iname(idi)

! -----

      end subroutine s_getiname

!-----7--------------------------------------------------------------7--

      end module m_getiname
