!***********************************************************************
      module m_upwnp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/04/03
!     Modification: 2006/09/30, 2007/10/19, 2007/11/26, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2009/11/05, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the sedimentation for optional precipitation
!     concentrations.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: upwnp, s_upwnp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface upwnp

        module procedure s_upwnp

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
      subroutine s_upwnp(fpdziv,dtp,ni,nj,nk,rbr,rst,un,ncf,ncflx)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtp
                       ! Time steps interval of fall out integration

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: un(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity
                       ! of optional precipitation concentrations

! Input and output variable

      real, intent(inout) :: ncf(0:ni+1,0:nj+1,1:nk)
                       ! Optional precipitation concentrations at future

! Internal shared variables

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

      real dziv        ! Inverse of dz

      real dzvdt       ! dziv x dtp

      real, intent(inout) :: ncflx(0:ni+1,0:nj+1,1:nk)
                       ! Fallout flux of
                       ! optional precipitation concentrations

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

      dzvdt=dziv*dtp

! -----

! Calculate the sedimentation.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ncflx(i,j,k)=rbr(i,j,k)*un(i,j,k)*ncf(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=1,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ncf(i,j,k)=max(ncf(i,j,k)                                     &
     &      +(ncflx(i,j,k+1)-ncflx(i,j,k))/rst(i,j,k)*dzvdt,0.e0)
        end do
        end do

!$omp end do

      end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        ncf(i,j,nkm1)=ncf(i,j,nkm2)
      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_upwnp

!-----7--------------------------------------------------------------7--

      end module m_upwnp
