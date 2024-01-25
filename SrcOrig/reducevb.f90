!***********************************************************************
      module m_reducevb
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/04/04
!     Modification: 2006/12/04, 2007/01/20, 2008/05/02, 2008/07/25,
!                   2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     reduce the integrated value on bottom and top boundary.

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

      public :: reducevb, s_reducevb

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface reducevb

        module procedure s_reducevb

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
      subroutine s_reducevb(sum0)
!***********************************************************************

! Input and output variable

      real, intent(inout) :: sum0
                       ! Integrated value on top and bottom boundary

! Internal shared variables

      integer ierr     ! Error descriptor

      real tmp1        ! Temporary variable

!-----7--------------------------------------------------------------7--

! Reduce the integrated value in each processor element and hold the
! common integrated value between all processor elements.

      call mpi_allreduce(sum0,tmp1,1,mpi_real,mpi_sum,mpi_comm_cress,   &
     &                   ierr)

      sum0=tmp1

! -----

      end subroutine s_reducevb

!-----7--------------------------------------------------------------7--

      end module m_reducevb
