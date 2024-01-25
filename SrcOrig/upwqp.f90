!***********************************************************************
      module m_upwqp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/04/03
!     Modification: 2006/05/12, 2006/09/30, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2009/11/05, 2011/03/18,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the sedimentation and precipitation for optional
!     precipitation mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: upwqp, s_upwqp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface upwqp

        module procedure s_upwqp

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
      subroutine s_upwqp(fpadvopt,fpdziv,dtp,ni,nj,nk,rbr,rst,uq,qpf,   &
     &                   precip,qpflx)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

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

      real, intent(in) :: uq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity
                       ! of optional precipitation mixing ratio

! Input and output variables

      real, intent(inout) :: qpf(0:ni+1,0:nj+1,1:nk)
                       ! Optional precipitation mixing ratio at future

      real, intent(inout) :: precip(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation
                       ! of optional precipitation mixing ratio

! Internal shared variables

      integer advopt   ! Option for advection scheme

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

      real dziv        ! Inverse of dz

      real dtp05       ! 0.5 x dtp

      real dzvdt       ! dziv x dtp

      real rwiv05      ! 0.5 / rhow

      real, intent(inout) :: qpflx(0:ni+1,0:nj+1,1:nk)
                       ! Fallout flux of
                       ! optional precipitation mixing ratio

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpadvopt,advopt)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

      dtp05=.5e0*dtp

      dzvdt=dziv*dtp

      rwiv05=.5e0/rhow

! -----

! Calculate the sedimentation and precipitation.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qpflx(i,j,k)=rbr(i,j,k)*uq(i,j,k)*qpf(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=1,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qpf(i,j,k)=max(qpf(i,j,k)                                     &
     &      +(qpflx(i,j,k+1)-qpflx(i,j,k))/rst(i,j,k)*dzvdt,0.e0)
        end do
        end do

!$omp end do

      end do

      if(advopt.le.3) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qpf(i,j,nkm1)=qpf(i,j,nkm2)

          precip(i,j,1)=(qpflx(i,j,1)+qpflx(i,j,2))*rwiv05
          precip(i,j,2)=precip(i,j,2)+precip(i,j,1)*dtp05

        end do
        end do

!$omp end do

      else

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qpf(i,j,nkm1)=qpf(i,j,nkm2)

          precip(i,j,1)=(qpflx(i,j,1)+qpflx(i,j,2))*rwiv05
          precip(i,j,2)=precip(i,j,2)+precip(i,j,1)*dtp

        end do
        end do

!$omp end do

      end if

!$omp end parallel

! -----

      end subroutine s_upwqp

!-----7--------------------------------------------------------------7--

      end module m_upwqp
