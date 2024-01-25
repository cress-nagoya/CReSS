!***********************************************************************
      module m_kh8uv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/06/10
!     Modification: 2006/02/13, 2006/11/06, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/08/09, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the horizontal eddy diffusivity at the u and v points.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: kh8uv, s_kh8uv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface kh8uv

        module procedure s_kh8uv

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
      subroutine s_kh8uv(fpmpopt,fpmfcopt,ni,nj,nk,rmf,rkh,rkh8u,rkh8v)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

! Input and output variable

      real, intent(inout) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy diffusivity / Jacobian

! Output variables

      real, intent(out) :: rkh8u(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity
                       ! / Jacobian at u points

      real, intent(out) :: rkh8v(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity
                       ! / Jacobian at v points

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

! Remark

!     rkh: This variable is also temporary, because it is not used
!          again.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)

! -----

! Set the horizontal eddy diffusivity at the u and v points.

!$omp parallel default(shared) private(k)

      if(mfcopt.eq.1.and.(mpopt.eq.0.or.mpopt.eq.10)) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-1
            rkh8u(i,j,k)=rkh(i-1,j,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            rkh(i,j,k)=rmf(i,j,2)*rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=1,ni-1
            rkh8v(i,j,k)=rkh(i,j-1,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

      else if(mfcopt.eq.1.and.mpopt.eq.5) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=1,ni-1
            rkh8v(i,j,k)=rkh(i,j-1,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            rkh(i,j,k)=rmf(i,j,2)*rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-1
            rkh8u(i,j,k)=rkh(i-1,j,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-1
            rkh8u(i,j,k)=rkh(i-1,j,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=1,ni-1
            rkh8v(i,j,k)=rkh(i,j-1,k)+rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! -----

      end subroutine s_kh8uv

!-----7--------------------------------------------------------------7--

      end module m_kh8uv
