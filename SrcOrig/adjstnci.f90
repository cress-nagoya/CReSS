!***********************************************************************
      module m_adjstnci
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2001/10/18, 2001/11/20, 2002/01/07,
!                   2002/01/15, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2003/12/12, 2004/04/01, 2004/05/31,
!                   2004/06/10, 2004/09/01, 2004/09/25, 2004/10/12,
!                   2004/12/17, 2005/04/04, 2005/09/30, 2005/10/05,
!                   2006/04/03, 2007/10/19, 2007/11/26, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2009/11/05, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     adjust the concentrations of the cloud ice.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: adjstnci, s_adjstnci

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface adjstnci

        module procedure s_adjstnci

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_adjstnci(ni,nj,nk,qi,nci)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

! Input and output variable

      real, intent(inout) :: nci(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice

! Internal shared variables

      real mimiv5      ! 0.5 / mimax
      real mi0iv2      ! 100.0 / mi0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      mimiv5=.5e0/mimax
      mi0iv2=1.e2/mi0

! -----

! Adjust the concentrations of the cloud ice.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          nci(i,j,k)                                                    &
     &      =min(max(nci(i,j,k),mimiv5*qi(i,j,k)),mi0iv2*qi(i,j,k))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_adjstnci

!-----7--------------------------------------------------------------7--

      end module m_adjstnci
