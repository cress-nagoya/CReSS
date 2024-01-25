!***********************************************************************
      module m_mk2cgbw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/09/30
!     Modification: 2007/03/29, 2007/10/19, 2008/01/11, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2009/02/27, 2011/08/18,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     change measurement from [m] [k] to [c] [g] for warm bin cloud
!     physics.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: mk2cgbw, s_mk2cgbw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface mk2cgbw

        module procedure s_mk2cgbw

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
      subroutine s_mk2cgbw(ni,nj,nk,nqw,nnw,rbr,rst,rbv,p,qwbin,nwbin,  &
     &                     prr)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

! Input and output variables

      real, intent(inout) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(inout) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(inout) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(inout) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(inout) :: qwbin(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water mixing ratio

      real, intent(inout) :: nwbin(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations

      real, intent(inout) :: prr(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for rain

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in bin categories

!-----7--------------------------------------------------------------7--

! Change measurement from [m] [k] to [c] [g].

!$omp parallel default(shared) private(k,n)

      do n=1,nqw

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            qwbin(i,j,k,n)=1.e-3*rbr(i,j,k)*qwbin(i,j,k,n)
            nwbin(i,j,k,n)=1.e-6*rbr(i,j,k)*nwbin(i,j,k,n)
          end do
          end do

!$omp end do

        end do

      end do

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          rbr(i,j,k)=1.e-3*rbr(i,j,k)
          rst(i,j,k)=1.e-3*rst(i,j,k)
          rbv(i,j,k)=1.e3*rbv(i,j,k)
          p(i,j,k)=1.e1*p(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        prr(i,j,1)=1.e2*prr(i,j,1)
        prr(i,j,2)=1.e2*prr(i,j,2)
      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_mk2cgbw

!-----7--------------------------------------------------------------7--

      end module m_mk2cgbw
