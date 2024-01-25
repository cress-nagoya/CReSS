!***********************************************************************
      module m_getunit
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/06/21, 2000/01/17,
!                   2001/04/15, 2001/05/29, 2002/06/18, 2003/04/30,
!                   2003/05/19, 2004/01/09, 2004/09/25, 2007/01/20,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the unit number to open the file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comionum

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getunit, s_getunit

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getunit

        module procedure s_getunit

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
      subroutine s_getunit(io)
!***********************************************************************

! Output variable

      integer, intent(out) :: io
                       ! i/o unit number

! Internal shared variable

      integer iio      ! Index of unit numbers table

!-----7--------------------------------------------------------------7--

! Get the unit number to open the file.

      io=iolst(1)

      do iio=1,nio-1
        iolst(iio)=iolst(iio+1)
      end do

      iolst(nio)=iolst(nio)+1

! -----

      end subroutine s_getunit

!-----7--------------------------------------------------------------7--

      end module m_getunit
