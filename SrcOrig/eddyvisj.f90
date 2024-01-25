!***********************************************************************
      module m_eddyvisj
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/02/13
!     Modification: 2006/11/06, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     the eddy viscosity is devided by Jacobian.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: eddyvisj, s_eddyvisj

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface eddyvisj

        module procedure s_eddyvisj

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
      subroutine s_eddyvisj(fpmfcopt,ni,nj,nk,jcb,mf,rkh,rkv)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

! Input and output variables

      real, intent(inout) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy diffusivity / Jacobian

      real, intent(inout) :: rkv(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy diffusivity / Jacobian

! Internal shared variable

      integer mfcopt   ! Option for map scale factor

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real jcbiv       ! Inverse of Jacobian

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpmfcopt,mfcopt)

! -----

! The eddy viscosity is devided by Jacobian.

!$omp parallel default(shared) private(k)

      if(mfcopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,jcbiv)

          do j=1,nj-1
          do i=1,ni-1
            jcbiv=1.e0/jcb(i,j,k)

            rkh(i,j,k)=jcbiv*rkh(i,j,k)
            rkv(i,j,k)=jcbiv*rkv(i,j,k)

          end do
          end do

!$omp end do

        end do

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,jcbiv)

          do j=1,nj-1
          do i=1,ni-1
            jcbiv=1.e0/jcb(i,j,k)

            rkh(i,j,k)=jcbiv*mf(i,j)*rkh(i,j,k)
            rkv(i,j,k)=jcbiv*rkv(i,j,k)

          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! -----

      end subroutine s_eddyvisj

!-----7--------------------------------------------------------------7--

      end module m_eddyvisj
