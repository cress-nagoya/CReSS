!***********************************************************************
      module m_convland
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/05/14
!     Modification: 2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     convert variable type of the land use categories.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: convland, s_convland

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface convland

        module procedure s_convland

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic nint
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_convland(fproc,ni,nj,land,rland)
!***********************************************************************

! Input variables

      character(len=7), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

! Input and output variables

      integer, intent(inout) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(inout) :: rland(0:ni+1,0:nj+1)
                       ! Real land use of surface

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

!! Convert variable type of the land use categories.

!$omp parallel default(shared)

! Convert the integer land use categories to real.

      if(fproc(1:4).eq.'real') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          rland(i,j)=real(land(i,j))+.1e0
        end do
        end do

!$omp end do

! -----

! Convert the real land use categories to integer.

      else if(fproc(1:7).eq.'integer') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          land(i,j)=nint(rland(i,j))
        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_convland

!-----7--------------------------------------------------------------7--

      end module m_convland
