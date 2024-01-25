!***********************************************************************
      module m_castvar
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2003/05/19, 2003/11/05, 2006/12/04, 2008/05/02,
!                   2008/07/25, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     broadcast optional variable from the root processor element to the
!     others.

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

      public :: castvar, s_castvar

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface castvar

        module procedure s_castvar

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
      subroutine s_castvar(imin,imax,jmin,jmax,kmin,kmax,var)
!***********************************************************************

! Input variable

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

! Input and output variables

      real, intent(inout) :: var(imin:imax,jmin:jmax,kmin:kmax)
                       ! Optional 3 dimensional variable

! Internal shared variables

      integer siz      ! Broadcasting buffer size

      integer ierr     ! Error descriptor

!-----7--------------------------------------------------------------7--

! Broadcast optional variable.

      siz=(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)

      call mpi_bcast(var,siz,mpi_real,root,mpi_comm_cress,ierr)

! -----

      end subroutine s_castvar

!-----7--------------------------------------------------------------7--

      end module m_castvar
