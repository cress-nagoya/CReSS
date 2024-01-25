!***********************************************************************
      module m_getrname
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/03/25
!     Modification: 2000/01/17, 2002/04/02, 2003/05/19, 2004/01/09,
!                   2007/01/20, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the real namelist variable.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getrname, s_getrname

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getrname

        module procedure s_getrname

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

      end module m_getrname
