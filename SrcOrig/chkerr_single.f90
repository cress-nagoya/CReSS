!***********************************************************************
      module m_chkerr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/06/18
!     Modification: 2003/05/19, 2008/05/02, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the error descriptor for calling procedure, destroy.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkerr, s_chkerr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkerr

        module procedure s_chkerr

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
      subroutine s_chkerr(stat)
!***********************************************************************

! Input and output variable

      integer, intent(inout) :: stat
                       ! Runtime status

!-----7--------------------------------------------------------------7--

! Set the error descriptor.

      if(stat.eq.0) then
        stat=mype+1
      else
        stat=-mype-1
      end if

! -----

      end subroutine s_chkerr

!-----7--------------------------------------------------------------7--

      end module m_chkerr
