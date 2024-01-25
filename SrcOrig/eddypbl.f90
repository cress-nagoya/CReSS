!***********************************************************************
      module m_eddypbl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/16
!     Modification: 2001/11/14, 2002/01/15, 2002/04/02, 2002/12/27,
!                   2003/02/13, 2003/03/13, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2004/02/01, 2004/03/05, 2004/04/01,
!                   2004/05/07, 2004/07/01, 2004/08/20, 2006/05/12,
!                   2007/06/27, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the eddy viscosity and diffusivity in planetaty boundary
!     layer.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: eddypbl, s_eddypbl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface eddypbl

        module procedure s_eddypbl

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_eddypbl(fplevpbl,fpdz,ni,nj,nk,zph,jcb8w,u,v,ptv,    &
     &                     kms,khs,vk)
!***********************************************************************

! Input variables

      integer, intent(in) :: fplevpbl
                       ! Formal parameter of unique index of levpbl

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature

! Output variables

      real, intent(out) :: kms(0:ni+1,0:nj+1,1:nk)
                       ! Eddy viscosity in planetaty boundary layer

      real, intent(out) :: khs(0:ni+1,0:nj+1,1:nk)
                       ! Eddy diffusivity in planetaty boundary layer

! Internal shared variables

      integer levpbl   ! Number of planetary boundary layer

      real dz          ! Grid distance in z direction

      real dzv225      ! 0.25 / (dz x dz)

      real gdziv2      ! 2.0 x g / dz

      real, intent(inout) :: vk(0:ni+1,0:nj+1,1:nk)
                       ! Square of virtical shear

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real grch        ! Gradient Richardson number
      real frch        ! Flux Richardson number

      real frchm1      ! 1.0 - frch

      real htwkp       ! kappa x terrain height from surface

      real ln          ! Turbulent length scale

      real sm          ! Temporary variable
      real sh          ! Temporary variable

      real a           ! Temporary variable
      real b           ! Temporary variable

! Remark

!     kms,khs: These variables are also temporary.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fplevpbl,levpbl)
      call getrname(fpdz,dz)

! -----

! Set the common used variables.

      dzv225=.25e0/(dz*dz)

      gdziv2=2.e0*g/dz

! -----

!!! Calculate the eddy viscosity and diffusivity in planetaty boundary
!!! layer.

!$omp parallel default(shared) private(k)

! Calculate the square of virtical shear.

      do k=2,levpbl+1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          kms(i,j,k)=u(i,j,k)+u(i+1,j,k)
          khs(i,j,k)=v(i,j,k)+v(i,j+1,k)
        end do
        end do

!$omp end do

      end do

      do k=3,levpbl+1

!$omp do schedule(runtime) private(i,j,a,b)

        do j=1,nj-1
        do i=1,ni-1
          a=kms(i,j,k)-kms(i,j,k-1)
          b=khs(i,j,k)-khs(i,j,k-1)

          vk(i,j,k)=max(vamin,                                          &
     &      ((a*a+b*b)/(jcb8w(i,j,k)*jcb8w(i,j,k)))*dzv225)

        end do
        end do

!$omp end do

      end do

! -----

!! Calculate the eddy viscosity and diffusivity.

      do k=3,levpbl+1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,grch,frch,frchm1,htwkp,ln,sm,sh,a)

        do j=1,nj-1
        do i=1,ni-1

! Calculate the gradient Richardson number.

          grch=max(gdziv2*(ptv(i,j,k)-ptv(i,j,k-1))                     &
     &      /(jcb8w(i,j,k)*vk(i,j,k)*(ptv(i,j,k-1)+ptv(i,j,k))),rchmin)

! -----

! Calculate the flux Richardson number.

          frch=.725e0*(.186e0+grch-sqrt(.0346e0+grch*grch-.316e0*grch))

          frchm1=1.e0-frch

! -----

! Calculate the factor sh and sm.

          sh=max(2.34e0*(.229333e0-1.07467e0*frch)/frchm1,0.e0)

          sm=sh*(.173333e0-.641333e0*frch)/(.229333e0-.918667e0*frch)

          a=sqrt(15.e0*frchm1*sm*vk(i,j,k))

          sh=a*sh
          sm=a*sm

! -----

! Calculate the turbulent length scale.

          htwkp=kappa*(zph(i,j,k)-zph(i,j,2))

          if(grch.lt.0.e0) then

            ln=100.e0*htwkp/(htwkp+100.e0)

          else

            ln=30.e0*htwkp/(htwkp+30.e0)

          end if

! -----

! Finally get the eddy viscosity and diffusivity.

          a=ln*ln

          kms(i,j,k)=a*sm
          khs(i,j,k)=a*sh

! -----

        end do
        end do

!$omp end do

      end do

!! -----

!$omp end parallel

!!! -----

      end subroutine s_eddypbl

!-----7--------------------------------------------------------------7--

      end module m_eddypbl
