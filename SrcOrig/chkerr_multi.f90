!***********************************************************************
      module m_chkerr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/06/18
!     Modification: 2003/02/13, 2003/05/19, 2003/11/05, 2004/01/09,
!                   2004/08/31, 2006/12/04, 2007/01/20, 2008/05/02,
!                   2008/07/25, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the error descriptor for calling procedure, destroy.

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

! Internal shared variables

      integer ierr     ! Error descriptor

      real rstat       ! Real runtime status

      real tmp1        ! Temporary variable

!-----7--------------------------------------------------------------7--

! Reduce the error descriptor in each processor element and hold the
! common error descriptor between all processor elements.

      if(stat.ne.0) then
        rstat=-1.e0
      else
        rstat=1.e0
      end if

      call mpi_allreduce(rstat,tmp1,1,mpi_real,mpi_min,mpi_comm_cress,  &
     &                   ierr)

      rstat=tmp1

! -----

! Set the common runtime status.

      if(rstat.gt.0.e0) then
        stat=1
      else
        stat=-1
      end if

! -----

      end subroutine s_chkerr

!-----7--------------------------------------------------------------7--

      end module m_chkerr
