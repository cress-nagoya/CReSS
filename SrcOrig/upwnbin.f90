!***********************************************************************
      module m_upwnbin
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/09/30
!     Modification: 2007/10/19, 2008/01/11, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/02/27, 2009/11/05, 2011/08/18,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the sedimentation for optional bin concentrations.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: upwnbin, s_upwnbin

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface upwnbin

        module procedure s_upwnbin

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
      subroutine s_upwnbin(fpdziv,ncp,dtp,ni,nj,nk,nn,rbr,rst,ubn,      &
     &                     nbin,nbflx)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: ncp
                       ! Index of current processed bin

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nn
                       ! Number of categories of concentrations

      real, intent(in) :: dtp
                       ! Time steps interval of fall out integration

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density [g/cm^3]

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian [g/cm^3]

      real, intent(in) :: ubn(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity
                       ! of optional bin concentrations [cm/s]

! Input and output variable

      real, intent(inout) :: nbin(0:ni+1,0:nj+1,1:nk,1:nn)
                       ! Optional bin concentrations [1/cm^3]

! Internal shared variables

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

      real dziv        ! Inverse of dz

      real dzvdt2      ! 0.01 x dziv x dtp

      real, intent(inout) :: nbflx(0:ni+1,0:nj+1,1:nk)
                       ! Fallout flux of optional bin concentrations

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

      dzvdt2=.01e0*dziv*dtp

! -----

! Calculate the sedimentation.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          nbflx(i,j,k)=rbr(i,j,k)*ubn(i,j,k)*nbin(i,j,k,ncp)
        end do
        end do

!$omp end do

      end do

      do k=1,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          nbin(i,j,k,ncp)=max(nbin(i,j,k,ncp)                           &
     &      +(nbflx(i,j,k+1)-nbflx(i,j,k))/rst(i,j,k)*dzvdt2,0.e0)
        end do
        end do

!$omp end do

      end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        nbin(i,j,nkm1,ncp)=nbin(i,j,nkm2,ncp)
      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_upwnbin

!-----7--------------------------------------------------------------7--

      end module m_upwnbin
