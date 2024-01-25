!***********************************************************************
      module m_sumbin
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/08/08
!     Modification: 2006/09/30, 2007/05/14, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/08/18, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the total mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: sumbin, s_sumbin

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface sumbin

        module procedure s_sumbin

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
      subroutine s_sumbin(ni,nj,nk,nq,rbv,mbin,qall)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nq
                       ! Number of categories of hydrometeor

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density [cm^3/g]

      real, intent(in) :: mbin(0:ni+1,0:nj+1,1:nk,1:nq)
                       ! Optional bin mass [g/cm^3]

! Output variable

      real, intent(out) :: qall(0:ni+1,0:nj+1,1:nk)
                       ! Total mixing ratio [g/g]

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in bin categories

!-----7--------------------------------------------------------------7--

! Get the total mixing ratio.

!$omp parallel default(shared) private(k,n)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qall(i,j,k)=rbv(i,j,k)*mbin(i,j,k,1)
        end do
        end do

!$omp end do

      end do

      do n=2,nq

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            qall(i,j,k)=qall(i,j,k)+rbv(i,j,k)*mbin(i,j,k,n)
          end do
          end do

!$omp end do

        end do

      end do

!$omp end parallel

! -----

      end subroutine s_sumbin

!-----7--------------------------------------------------------------7--

      end module m_sumbin
