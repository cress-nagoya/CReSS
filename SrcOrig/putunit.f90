!***********************************************************************
      module m_putunit
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 2000/01/17, 2001/04/15,
!                   2003/04/30, 2003/05/19, 2004/01/09, 2004/09/25,
!                   2007/01/20, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2009/11/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     return the unit number to close the file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comionum

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: putunit, s_putunit

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface putunit

        module procedure s_putunit

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
      subroutine s_putunit(io)
!***********************************************************************

! Input variable

      integer, intent(in) :: io
                       ! i/o unit number

! Internal shared variables

      integer iio      ! Index of unit numbers table

      integer insert   ! Returned unit number

!-----7--------------------------------------------------------------7--

! Return the unit number to close the file.

      if(io.lt.iolst(1)) then

        insert=1

      else

        insert=nio

        do iio=1,nio-1

          if(io.gt.iolst(iio)) then

            if(iio.lt.insert) then
              insert=iio
            end if

          end if

        end do

        insert=insert+1

      end if

      do iio=nio,insert+1,-1
        iolst(iio)=iolst(iio-1)
      end do

      iolst(insert)=io

! -----

      end subroutine s_putunit

!-----7--------------------------------------------------------------7--

      end module m_putunit
