!***********************************************************************
      module m_uvw2gpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/03/13
!     Modification: 2001/04/15, 2001/06/29, 2002/04/02, 2002/08/15,
!                   2002/09/09, 2002/10/31, 2003/04/30, 2003/05/19,
!                   2003/11/28, 2003/12/12, 2006/09/21, 2007/05/07,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/03/23, 2011/09/22, 2013/01/28, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the analysis nudging to GPV data of the velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: uvw2gpv, s_uvw2gpv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface uvw2gpv

        module procedure s_uvw2gpv

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
      subroutine s_uvw2gpv(fpnggvar,nggdmp,gtinc,                       &
     &                     ni,nj,nk,rst8u,rst8v,rst8w,up,vp,wp,         &
     &                     ugpv,utd,vgpv,vtd,wgpv,wtd,ufrc,vfrc,wfrc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpnggvar
                       ! Formal parameter of unique index of nggvar

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) ::  rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jabobian at u points

      real, intent(in) ::  rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jabobian at v points

      real, intent(in) ::  rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jabobian at w points

      real, intent(in) ::  up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) ::  vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) ::  wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) ::  ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(in) ::  utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(in) ::  vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(in) ::  vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(in) ::  wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(in) ::  wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

! Input and output variables

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variable

      character(len=108) nggvar
                       ! Control flag of
                       ! analysis nudged variables to GPV

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(nggvar)

! -----

! Get the required namelist variable.

      call getcname(fpnggvar,nggvar)

! -----

!! Perform the analysis nudging to GPV data of the velocity.

!$omp parallel default(shared) private(k)

! Calculate the analysis nudging term for the x components of velocity.

      if(nggvar(1:1).eq.'o') then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            ufrc(i,j,k)=ufrc(i,j,k)+nggdmp                              &
     &        *rst8u(i,j,k)*((ugpv(i,j,k)+utd(i,j,k)*gtinc)-up(i,j,k))
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the analysis nudging term for the y components of velocity.

      if(nggvar(2:2).eq.'o') then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            vfrc(i,j,k)=vfrc(i,j,k)+nggdmp                              &
     &        *rst8v(i,j,k)*((vgpv(i,j,k)+vtd(i,j,k)*gtinc)-vp(i,j,k))
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the analysis nudging term for the z components of velocity.

      if(nggvar(3:3).eq.'o') then

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wfrc(i,j,k)=wfrc(i,j,k)+nggdmp                              &
     &        *rst8w(i,j,k)*((wgpv(i,j,k)+wtd(i,j,k)*gtinc)-wp(i,j,k))
          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_uvw2gpv

!-----7--------------------------------------------------------------7--

      end module m_uvw2gpv
