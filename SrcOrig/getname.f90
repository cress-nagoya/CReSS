!***********************************************************************
      module m_getname
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2013/02/13
!     Modification: 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the namelist variable.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getname, s_getcname, s_getiname, s_getrname

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getname

        module procedure s_getcname, s_getiname, s_getrname

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

!***********************************************************************
      subroutine s_getrname(idr,rval)
!***********************************************************************

! Input variable

      integer, intent(in) :: idr
                       ! Optional index in real namelist table

! Output variable

      real, intent(out) :: rval
                       ! Real namelist variable

!-----7--------------------------------------------------------------7--

! Get the real namelist variable.

      rval=rname(idr)

! -----

      end subroutine s_getrname

!-----7--------------------------------------------------------------7--

      end module m_getname
