!***********************************************************************
      module m_chkstd
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/06/18
!     Modification: 2003/05/19, 2004/05/31, 2004/08/20, 2006/12/04,
!                   2007/01/20, 2008/05/02, 2008/07/25, 2008/08/25,
!                   2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the control flag of message type for standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_comstd
      use m_defmpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkstd, s_chkstd

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkstd

        module procedure s_chkstd

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
      subroutine s_chkstd(broot)
!***********************************************************************

! Input variable

      integer, intent(in) :: broot
                       ! Broadcasting root

! Internal shared variable

      integer ierr     ! Error descriptor

!-----7--------------------------------------------------------------7--

! Hold the common value of the control flag.

      call mpi_bcast(fstd,3,mpi_character,broot,mpi_comm_cress,ierr)

! -----

      end subroutine s_chkstd

!-----7--------------------------------------------------------------7--

      end module m_chkstd
