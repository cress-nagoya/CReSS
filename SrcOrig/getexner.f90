!***********************************************************************
      module m_getexner
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/03/08
!     Modification: 2000/04/18, 2000/07/05, 2001/05/29, 2002/04/02,
!                   2003/04/30, 2003/05/19, 2004/04/15, 2007/10/19,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the total pressure variable and Exner function.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getexner, s_getexner

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getexner

        module procedure s_getexner

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getexner(ni,nj,nk,pbr,pp,pi,p)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

! Output variables

      real, intent(out) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(out) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exner function

! Internal shared variables

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      rddvcp=rd/cp

      p0iv=1.e0/p0

! -----

! Calculate the total pressure variable and Exner function.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          p(i,j,k)=pbr(i,j,k)+pp(i,j,k)

          pi(i,j,k)=exp(rddvcp*log(p0iv*p(i,j,k)))

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_getexner

!-----7--------------------------------------------------------------7--

      end module m_getexner
