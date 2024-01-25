!***********************************************************************
      module m_chkitr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/06/01
!     Modification: 2000/07/05, 2001/04/15, 2001/05/29, 2001/08/07,
!                   2002/04/02, 2003/04/30, 2003/05/19, 2003/11/05,
!                   2004/08/20, 2004/09/10, 2006/12/04, 2007/01/20,
!                   2007/05/14, 2007/10/19, 2008/05/02, 2008/07/25,
!                   2008/08/25, 2009/02/27, 2009/11/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the convergence of the iteration.

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

      public :: chkitr, s_chkitr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkitr

        module procedure s_chkitr

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chkitr(fproc,istr,iend,jstr,jend,kstr,kend,itcon,    &
     &                    ni,nj,nk,dvar)
!***********************************************************************

! Input variables

      character(len=6), intent(in) :: fproc
                       ! Control flag of parallel processing

      integer, intent(in) :: istr
                       ! Minimum do loops index in x direction

      integer, intent(in) :: iend
                       ! Maximum do loops index in x direction

      integer, intent(in) :: jstr
                       ! Minimum do loops index in y direction

      integer, intent(in) :: jend
                       ! Maximum do loops index in y direction

      integer, intent(in) :: kstr
                       ! Minimum do loops index in z direction

      integer, intent(in) :: kend
                       ! Maximum do loops index in z direction

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dvar(0:ni+1,0:nj+1,1:nk)
                       ! Variations of iteration

! Output variable

      real, intent(out) :: itcon
                       ! Control flag of continuation of interation

! Internal shared variables

      integer ierr     ! Error descriptor

      real intitc      ! Internal processed control flag of
                       ! continuation of interation

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Initialize the processed variable, intitc.

      intitc=0.e0

! -----

! Check the convergence of the iteration in each processor element.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,k) reduction(max: intitc)

      do k=kstr,kend
      do j=jstr,jend
      do i=istr,iend
        intitc=max(abs(dvar(i,j,k)),intitc)
      end do
      end do
      end do

!$omp end do

!$omp end parallel

! -----

! Hold the common control flag itcon.

      if(fproc(1:6).eq.'common') then

        call mpi_allreduce(intitc,itcon,1,mpi_real,mpi_max,             &
     &                     mpi_comm_cress,ierr)

      end if

! -----

      end subroutine s_chkitr

!-----7--------------------------------------------------------------7--

      end module m_chkitr
