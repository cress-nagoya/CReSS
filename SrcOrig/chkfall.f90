!***********************************************************************
      module m_chkfall
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2001/06/29, 2003/05/19, 2003/11/05, 2004/09/25,
!                   2006/04/03, 2007/01/20, 2008/05/02, 2008/07/25,
!                   2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the time interval of upwind differential of sedimentation to
!     hold the common value.

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

      public :: chkfall, s_chkfall

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkfall

        module procedure s_chkfall

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
      subroutine s_chkfall(delt)
!***********************************************************************

! Input and output variable

      real, intent(inout) :: delt
                       ! Time interval of upwind differential

! Internal shared variables

      integer ierr     ! Error descriptor

      real tmp1        ! Temporary variable

!-----7--------------------------------------------------------------7--

! Reduce the time interval in each processor element and hold the common
! time interval between all processor elements.

      call mpi_allreduce(delt,tmp1,1,mpi_real,mpi_min,mpi_comm_cress,   &
     &                   ierr)

      delt=tmp1

! -----

      end subroutine s_chkfall

!-----7--------------------------------------------------------------7--

      end module m_chkfall
