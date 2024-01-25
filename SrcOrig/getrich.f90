!***********************************************************************
      module m_getrich
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/15
!     Modification: 2001/10/18, 2001/11/14, 2002/01/15, 2002/04/02,
!                   2002/07/03, 2002/12/02, 2002/12/27, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/07/15, 2003/10/31,
!                   2003/12/12, 2004/03/05, 2004/04/01, 2004/05/07,
!                   2004/08/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/11/13, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the bulk Richardson number on the surface.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getrich, s_getrich

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getrich

        module procedure s_getrich

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
      subroutine s_getrich(ni,nj,nk,za,land,kai,z0m,z0h,ptv,va,rch)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: za(0:ni+1,0:nj+1)
                       ! z physical coordinates at lowest plane

      real, intent(in) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

      real, intent(in) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(in) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(in) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature

      real, intent(in) :: va(0:ni+1,0:nj+1)
                       ! Square root of velocity at lowest plane

! Output variable

      real, intent(out) :: rch(0:ni+1,0:nj+1)
                       ! Bulk Richardson number

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      real dz0m        ! za - z0m
      real dz0h        ! za - z0h

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Calculate the bulk Richardson number.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j,dz0m,dz0h,a)

      do j=1,nj-1
      do i=1,ni-1
        dz0m=za(i,j)-z0m(i,j)
        dz0h=za(i,j)-z0h(i,j)

        a=g*(ptv(i,j,2)-ptv(i,j,1))/(ptv(i,j,1)*va(i,j)*va(i,j))

        rch(i,j)=max(a*dz0m*dz0m/dz0h,rchmin)

        if(land(i,j).eq.1) then

          dz0m=za(i,j)-icz0m
          dz0h=za(i,j)-icz0h

          rch(i,j)=(1.e0-kai(i,j))*rch(i,j)                             &
     &      +kai(i,j)*max(a*dz0m*dz0m/dz0h,rchmin)

        end if

      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_getrich

!-----7--------------------------------------------------------------7--

      end module m_getrich
