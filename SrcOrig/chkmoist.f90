!***********************************************************************
      module m_chkmoist
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 2000/01/17, 2001/04/15, 2001/05/29, 2002/04/02,
!                   2002/06/06, 2002/08/15, 2003/04/30, 2003/05/19,
!                   2003/11/05, 2004/01/09, 2004/08/31, 2004/09/10,
!                   2006/12/04, 2007/01/20, 2007/05/07, 2007/05/14,
!                   2007/10/19, 2008/05/02, 2008/07/25, 2008/08/25,
!                   2009/02/27, 2009/11/05

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the air moisture.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_commpi
      use m_defmpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkmoist, s_chkmoist

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkmoist

        module procedure s_chkmoist

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chkmoist(fmois,ni,nj,nk,qv)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

! Output variable

      character(len=5), intent(out) :: fmois
                       ! Control flag of air moisture

! Internal shared variables

      integer ierr     ! Error descriptor

      real qvmax       ! Maximum water vapor mixing ratio

      real tmp1        ! Temporary variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Initialize the processed variable, qvmax.

      qvmax=lim36n

! -----

! Check the air moisture in each processor element.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,k) reduction(max: qvmax)

      do k=1,nk-1
      do j=1,nj-1
      do i=1,ni-1
        qvmax=max(qv(i,j,k),qvmax)
      end do
      end do
      end do

!$omp end do

!$omp end parallel

! -----

! Reduce the maximum value in each processor element and hold the common
! maximum value between all processor elements.

      call mpi_allreduce(qvmax,tmp1,1,mpi_real,mpi_max,mpi_comm_cress,  &
     &                   ierr)

      qvmax=tmp1

! -----

! Set the common control flag fmois.

      if(qvmax.gt.0.e0) then

        write(fmois(1:5),'(a5)') 'moist'

      else

        write(fmois(1:5),'(a5)') 'dry  '

      end if

! -----

      end subroutine s_chkmoist

!-----7--------------------------------------------------------------7--

      end module m_chkmoist
