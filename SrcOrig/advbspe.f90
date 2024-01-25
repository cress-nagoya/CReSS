!***********************************************************************
      module m_advbspe
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/08/03, 1999/09/30, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/07/05, 2002/04/02, 2003/04/30,
!                   2003/05/19, 2003/12/26, 2005/08/05, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the base state pressure advection for the horizontally
!     explicit and vertically explicit method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: advbspe, s_advbspe

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface advbspe

        module procedure s_advbspe

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
      subroutine s_advbspe(ni,nj,nk,rst,w,pdiv,psml)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure x Jacobian

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: pdiv(0:ni+1,0:nj+1,1:nk)
                       ! Divergence value in pressure equation

! Output variable

      real, intent(out) :: psml(0:ni+1,0:nj+1,1:nk)
                       ! Acoustic mood in pressure equation

! Internal shared variable

      real g05         ! 0.5 x g

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Set the common used variable.

      g05=.5e0*g

! -----

! Calculate the base state pressure advection.

!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          psml(i,j,k)=g05*(w(i,j,k)+w(i,j,k+1))*rst(i,j,k)+pdiv(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_advbspe

!-----7--------------------------------------------------------------7--

      end module m_advbspe
