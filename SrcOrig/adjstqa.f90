!***********************************************************************
      module m_adjstqa
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     force the aerosol mixing ratio more than user specified value.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: adjstqa, s_adjstqa

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface adjstqa

        module procedure s_adjstqa

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
      subroutine s_adjstqa(ni,nj,nk,nqa,qasl)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

! Input and output variable

      real, intent(inout) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in 4th direction

!-----7--------------------------------------------------------------7--

! Force the aerosol mixing ratio more than user specified value.

!$omp parallel default(shared) private(k,n)

      do n=1,nqa(0)

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            qasl(i,j,k,n)=max(qasl(i,j,k,n),0.e0)
          end do
          end do

!$omp end do

        end do

      end do

!$omp end parallel

! -----

      end subroutine s_adjstqa

!-----7--------------------------------------------------------------7--

      end module m_adjstqa
