!***********************************************************************
      module m_reducelb
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/10/28
!     Modification: 2003/05/19, 2003/11/05, 2006/12/04, 2007/01/20,
!                   2008/05/02, 2008/07/25, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     reduce the integrated value on lateral boundary.

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

      public :: reducelb, s_reducelb

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface reducelb

        module procedure s_reducelb

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
      subroutine s_reducelb(sumw,sume,sums,sumn)
!***********************************************************************

! Input and output variables

      real, intent(inout) :: sumw
                       ! Integrated value on west boundary

      real, intent(inout) :: sume
                       ! Integrated value on east boundary

      real, intent(inout) :: sums
                       ! Integrated value on south boundary

      real, intent(inout) :: sumn
                       ! Integrated value on north boundary

! Internal shared variables

      integer ierr     ! Error descriptor

      real ssum(1:4)   ! Sending buffer of sumw, sume, sums and sumn
      real rsum(1:4)   ! Receiving buffer of sumw, sume, sums and sumn

!-----7--------------------------------------------------------------7--

! Put the integrated value in each processor element to the sending
! buffer.

      ssum(1)=sumw
      ssum(2)=sume
      ssum(3)=sums
      ssum(4)=sumn

! -----

! Reduce the integrated value in each processor element and hold the
! common integrated value between all processor elements.

      call mpi_allreduce(ssum,rsum,4,mpi_real,mpi_sum,mpi_comm_cress,   &
     &                   ierr)

! -----

! Get the common integrated value from the receiving buffer.

      sumw=rsum(1)
      sume=rsum(2)
      sums=rsum(3)
      sumn=rsum(4)

! -----

      end subroutine s_reducelb

!-----7--------------------------------------------------------------7--

      end module m_reducelb
