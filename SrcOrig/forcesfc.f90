!***********************************************************************
      module m_forcesfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/03/13
!     Modification: 2003/04/30, 2003/05/19, 2004/02/01, 2004/03/05,
!                   2004/09/10, 2005/06/10, 2006/01/10, 2007/01/20,
!                   2007/05/21, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/09/22, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the surface flux to bottom boundary.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forcesfc, s_forcesfc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forcesfc

        module procedure s_forcesfc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_forcesfc(fmois,ni,nj,nk,j31,j32,ptbr,u,v,w,ptp,qv,   &
     &                      ptv,qvsfc,ce,ct,cq,ufrc,vfrc,ptfrc,qvfrc)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature

      real, intent(in) :: qvsfc(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio on surface

      real, intent(in) :: ce(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface momentum flux

      real, intent(in) :: ct(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface heat flux

      real, intent(in) :: cq(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface moisture flux

! Output variables

      real, intent(out) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(out) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(out) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in potential temperature equation

      real, intent(out) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real j318u       ! j31 at u point
      real j328v       ! j32 at v point

      real xcomp       ! Temporary variable
      real ycomp       ! Temporary variable
      real zcomp       ! Temporary variable

!-----7--------------------------------------------------------------7--

!! Get the surface flux to bottom boundary.

!$omp parallel default(shared)

! Get the surface flux for the potential tempeture.

      if(fmois(1:3).eq.'dry') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ptfrc(i,j,1)=ct(i,j)*(ptv(i,j,2)-ptv(i,j,1))
        end do
        end do

!$omp end do

      else if(fmois(1:5).eq.'moist') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ptfrc(i,j,1)=ct(i,j)*((ptbr(i,j,2)+ptp(i,j,2))                &
     &      -ptv(i,j,1)*(1.e0+qvsfc(i,j))/(1.e0+epsav*qvsfc(i,j)))
        end do
        end do

!$omp end do

      end if

! -----

! Get the surface flux for water vapor mixing ratio.

      if(fmois(1:3).eq.'dry') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qvfrc(i,j,1)=0.e0
        end do
        end do

!$omp end do

      else if(fmois(1:5).eq.'moist') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qvfrc(i,j,1)=cq(i,j)*(qv(i,j,2)-qvsfc(i,j))
        end do
        end do

!$omp end do

      end if

! -----

! Get the surface flux for the x components of velocity.

!$omp do schedule(runtime) private(i,j,j318u,xcomp,zcomp)

      do j=1,nj-1
      do i=2,ni-1
        j318u=j31(i,j,2)+j31(i,j,3)

        xcomp=1.e0/sqrt(4.e0+j318u*j318u)
        zcomp=.125e0*j318u*xcomp

        ufrc(i,j,1)=(ce(i-1,j)+ce(i,j))*(u(i,j,2)*xcomp                 &
     &    +((w(i-1,j,2)+w(i,j,3))+(w(i-1,j,3)+w(i,j,2)))*zcomp)

      end do
      end do

!$omp end do

! -----

! Get the surface flux for the x components of velocity.

!$omp do schedule(runtime) private(i,j,j328v,ycomp,zcomp)

      do j=2,nj-1
      do i=1,ni-1
        j328v=j32(i,j,2)+j32(i,j,3)

        ycomp=1.e0/sqrt(4.e0+j328v*j328v)
        zcomp=.125e0*j328v*ycomp

        vfrc(i,j,1)=(ce(i,j-1)+ce(i,j))*(v(i,j,2)*ycomp                 &
     &    +((w(i,j-1,2)+w(i,j,3))+(w(i,j-1,3)+w(i,j,2)))*zcomp)

      end do
      end do

!$omp end do

! -----

!$omp end parallel

!! -----

      end subroutine s_forcesfc

!-----7--------------------------------------------------------------7--

      end module m_forcesfc
