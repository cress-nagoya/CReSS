!***********************************************************************
      module m_numchar
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 2000/01/17, 2001/04/15, 2003/04/30, 2003/05/19,
!                   2006/09/21, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     count the number of the optional character variable.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: numchar, s_numchar

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface numchar

        module procedure s_numchar

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
      subroutine s_numchar(opchar,sic,nc,ncspc)
!***********************************************************************

! Input variables

      character(len=108), intent(in) :: opchar
                       ! Optional character variable

      integer, intent(in) :: sic
                       ! Start point to count

! Output variables

      integer, intent(out) :: nc
                       ! Number of optional character

      integer, intent(out) :: ncspc
                       ! Number of space of optional character

! Internal shared variables

      integer ic       ! Index of optional character variable
      integer jc       ! Index of optional character variable

      integer ido      ! Control flag of continuation of increment

!-----7--------------------------------------------------------------7--

! Count the number of the optional character variable.

      ido=1

      ncspc=0

      do ic=sic,108

        if(ido.eq.1) then

          if(opchar(ic:ic).ne.' ') then
            ido=0
          else
            ncspc=ncspc+1
          end if

        end if

      end do

      ido=1

      nc=0

      do jc=ncspc+1,108

        if(ido.eq.1) then

          if(opchar(jc:jc).eq.' ') then
            ido=0
          else
            nc=nc+1
          end if

        end if

      end do

! -----

      end subroutine s_numchar

!-----7--------------------------------------------------------------7--

      end module m_numchar
