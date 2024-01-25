!***********************************************************************
      module m_coriuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/04/06, 1999/07/05, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/04/18, 2001/11/20, 2002/04/02,
!                   2002/11/11, 2003/01/04, 2003/04/30, 2003/05/19,
!                   2006/11/06, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the Coriolis force in the x, the y and the z components
!     of velocity equation.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: coriuvw, s_coriuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface coriuvw

        module procedure s_coriuvw

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
      subroutine s_coriuvw(ni,nj,nk,fc,rst,u,v,w,ufrc,vfrc,wfrc,        &
     &                     tmp1,tmp2)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: fc(0:ni+1,0:nj+1,1:2)
                       ! 0.25 x Coriolis parameters

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

! Input and output variables

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

!!! Calculate the Coriolis force in the x, the y and the z components
!!! of velocity equation.

!$omp parallel default(shared) private(k)

! Calculate the Coriolis force in the x components of velocity equation.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=1,ni-1
          tmp1(i,j,k)=rst(i,j,k)*(fc(i,j,1)*(v(i,j,k)+v(i,j+1,k))       &
     &      -fc(i,j,2)*(w(i,j,k)+w(i,j,k+1)))
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          ufrc(i,j,k)=ufrc(i,j,k)+(tmp1(i-1,j,k)+tmp1(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

!! Calculate the Coriolis force in the y and the z components of
!! velocity equation.

! Set the common used variable.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=2,ni-2
          tmp1(i,j,k)=rst(i,j,k)*(u(i,j,k)+u(i+1,j,k))
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the Coriolis force in the y components of velocity equation.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=2,ni-2
          tmp2(i,j,k)=-fc(i,j,1)*tmp1(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          vfrc(i,j,k)=vfrc(i,j,k)+(tmp2(i,j-1,k)+tmp2(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the Coriolis force in the z components of velocity equation.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          tmp2(i,j,k)=fc(i,j,2)*tmp1(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          wfrc(i,j,k)=wfrc(i,j,k)+(tmp2(i,j,k-1)+tmp2(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

!! -----

!$omp end parallel

!!! -----

      end subroutine s_coriuvw

!-----7--------------------------------------------------------------7--

      end module m_coriuvw
