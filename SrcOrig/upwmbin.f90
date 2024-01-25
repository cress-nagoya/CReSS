!***********************************************************************
      module m_upwmbin
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/09/30
!     Modification: 2007/10/19, 2008/01/11, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/02/27, 2009/11/05, 2011/03/18,
!                   2011/08/18, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the sedimentation and precipitation for optional bin
!     mass.

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

      public :: upwmbin, s_upwmbin

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface upwmbin

        module procedure s_upwmbin

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
      subroutine s_upwmbin(fpadvopt,fpdziv,ncp,dtp,ni,nj,nk,nq,         &
     &                     rbr,rst,rbv,ubm,mbin,precip,mbflx)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

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

      integer, intent(in) :: nq
                       ! Number of categories of hydrometeor

      real, intent(in) :: dtp
                       ! Time steps interval of fall out integration

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density [g/cm^3]

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian [g/cm^3]

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density [cm^3/g]

      real, intent(in) :: ubm(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of optional bin mass [cm/s]

! Input and output variables

      real, intent(inout) :: mbin(0:ni+1,0:nj+1,1:nk,1:nq)
                       ! Optional bin mass [g/cm^3]

      real, intent(inout) :: precip(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation of
                       ! optional bin mass [cm/s], [cm]

! Internal shared variables

      integer advopt   ! Option for advection scheme

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

      real dziv        ! Inverse of dz

      real dtp05       ! 0.5 x dtp

      real dzvdt2      ! 0.01 x dziv x dtp

      real rwv500      ! 500.0 / rhow

      real, intent(inout) :: mbflx(0:ni+1,0:nj+1,1:nk)
                       ! Fallout flux of optional precipitation mass

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

      dzvdt2=.01e0*dziv*dtp

      rwv500=500.e0/rhow

! -----

! Calculate the sedimentation and precipitation.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          mbflx(i,j,k)=rbr(i,j,k)*ubm(i,j,k)*mbin(i,j,k,ncp)
        end do
        end do

!$omp end do

      end do

      do k=1,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          mbin(i,j,k,ncp)=max(mbin(i,j,k,ncp)                           &
     &      +(mbflx(i,j,k+1)-mbflx(i,j,k))/rst(i,j,k)*dzvdt2,0.e0)
        end do
        end do

!$omp end do

      end do

      if(advopt.le.3) then

        if(ncp.eq.1) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            mbin(i,j,nkm1,ncp)=mbin(i,j,nkm2,ncp)

            precip(i,j,1)=rbv(i,j,1)*(mbflx(i,j,1)+mbflx(i,j,2))*rwv500

          end do
          end do

!$omp end do

        else if(ncp.eq.nq) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            mbin(i,j,nkm1,ncp)=mbin(i,j,nkm2,ncp)

            precip(i,j,1)=precip(i,j,1)                                 &
     &        +rbv(i,j,1)*(mbflx(i,j,1)+mbflx(i,j,2))*rwv500

            precip(i,j,2)=precip(i,j,2)+precip(i,j,1)*dtp05

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            mbin(i,j,nkm1,ncp)=mbin(i,j,nkm2,ncp)

            precip(i,j,1)=precip(i,j,1)                                 &
     &        +rbv(i,j,1)*(mbflx(i,j,1)+mbflx(i,j,2))*rwv500

          end do
          end do

!$omp end do

        end if

      else

        if(ncp.eq.1) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            mbin(i,j,nkm1,ncp)=mbin(i,j,nkm2,ncp)

            precip(i,j,1)=rbv(i,j,1)*(mbflx(i,j,1)+mbflx(i,j,2))*rwv500

          end do
          end do

!$omp end do

        else if(ncp.eq.nq) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            mbin(i,j,nkm1,ncp)=mbin(i,j,nkm2,ncp)

            precip(i,j,1)=precip(i,j,1)                                 &
     &        +rbv(i,j,1)*(mbflx(i,j,1)+mbflx(i,j,2))*rwv500

            precip(i,j,2)=precip(i,j,2)+precip(i,j,1)*dtp

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            mbin(i,j,nkm1,ncp)=mbin(i,j,nkm2,ncp)

            precip(i,j,1)=precip(i,j,1)                                 &
     &        +rbv(i,j,1)*(mbflx(i,j,1)+mbflx(i,j,2))*rwv500

          end do
          end do

!$omp end do

        end if

      end if

!$omp end parallel

! -----

      end subroutine s_upwmbin

!-----7--------------------------------------------------------------7--

      end module m_upwmbin
