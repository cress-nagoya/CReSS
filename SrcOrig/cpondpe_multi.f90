!***********************************************************************
      module m_cpondpe
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/25
!     Modification: 2000/01/17, 2001/02/13, 2001/04/15, 2002/07/03,
!                   2003/05/19, 2006/12/04, 2007/01/20, 2008/05/02,
!                   2008/07/25, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     wait for the other processor elements' signals.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_defmpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: cpondpe, s_cpondpe

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface cpondpe

        module procedure s_cpondpe

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
      subroutine s_cpondpe
!***********************************************************************

! Internal shared variable

      integer ierr     ! Error descriptor

!-----7--------------------------------------------------------------7--

! Wait for the other processor elements' signals.

      call mpi_barrier(mpi_comm_cress,ierr)

! -----

      end subroutine s_cpondpe

!-----7--------------------------------------------------------------7--

      end module m_cpondpe
