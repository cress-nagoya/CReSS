!***********************************************************************
      module m_vspuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/07/05, 1999/08/03, 1999/09/14, 1999/09/30,
!                   1999/10/12, 1999/11/01, 2000/01/17, 2001/06/29,
!                   2001/11/20, 2002/04/02, 2002/08/15, 2003/04/30,
!                   2003/05/19, 2003/12/12, 2004/04/15, 2006/01/10,
!                   2006/09/21, 2006/11/27, 2007/05/07, 2007/07/30,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2009/03/23, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the vertical sponge damping for the velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_getiname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vspuvw, s_vspuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vspuvw

        module procedure s_vspuvw

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
      subroutine s_vspuvw(fpgpvvar,fpvspvar,fpvspopt,ksp0,gtinc,        &
     &                    ni,nj,nk,ubr,vbr,rst8u,rst8v,rst8w,up,vp,wp,  &
     &                    rbct,ugpv,utd,vgpv,vtd,wgpv,wtd,              &
     &                    ufrc,vfrc,wfrc,rbct8s)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpvspvar
                       ! Formal parameter of unique index of vspvar

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(in) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

      real, intent(in) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(in) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(in) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

! Input and output variables

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) vspvar
                       ! Control flag of
                       ! vertical sponge damped variables

      integer vspopt   ! Option for vertical sponge damping

      real, intent(inout) :: rbct8s(1:ni,1:nj,1:nk)
                       ! 0.5 x relaxed top sponge damping coefficients
                       ! at scalar points

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(gpvvar)
      call inichar(vspvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getcname(fpvspvar,vspvar)
      call getiname(fpvspopt,vspopt)

! -----

!!! Calculate the vertical sponge damping for the velocity.

!$omp parallel default(shared) private(k)

! Set the common used variable.

      if(vspvar(1:1).eq.'o'.or.vspvar(2:2).eq.'o') then

        if(vspopt.eq.1) then

          do k=ksp0(1)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              rbct8s(i,j,k)=.25e0*(rbct(i,j,k,1)+rbct(i,j,k+1,1))
            end do
            end do

!$omp end do

          end do

        else

          do k=ksp0(2)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              rbct8s(i,j,k)=.25e0*(rbct(i,j,k,2)+rbct(i,j,k+1,2))
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!! For the x components of velocity.

      if(vspvar(1:1).eq.'o') then

! Damp to the GPV data.

        if(vspopt.eq.1) then

          do k=ksp0(1)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              ufrc(i,j,k)=ufrc(i,j,k)                                   &
     &          -rst8u(i,j,k)*(rbct8s(i-1,j,k)+rbct8s(i,j,k))           &
     &          *(up(i,j,k)-(ugpv(i,j,k)+utd(i,j,k)*gtinc))
            end do
            end do

!$omp end do

          end do

! -----

! Damp to the base state value.

        else

          do k=ksp0(2)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              ufrc(i,j,k)=ufrc(i,j,k)-rst8u(i,j,k)                      &
     &          *(rbct8s(i-1,j,k)+rbct8s(i,j,k))*(up(i,j,k)-ubr(i,j,k))
            end do
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!! For the y components of velocity.

      if(vspvar(2:2).eq.'o') then

! Damp to the GPV data.

        if(vspopt.eq.1) then

          do k=ksp0(1)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vfrc(i,j,k)=vfrc(i,j,k)                                   &
     &          -rst8v(i,j,k)*(rbct8s(i,j-1,k)+rbct8s(i,j,k))           &
     &          *(vp(i,j,k)-(vgpv(i,j,k)+vtd(i,j,k)*gtinc))
            end do
            end do

!$omp end do

          end do

! -----

! Damp to the base state value.

        else

          do k=ksp0(2)-1,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vfrc(i,j,k)=vfrc(i,j,k)-rst8v(i,j,k)                      &
     &          *(rbct8s(i,j-1,k)+rbct8s(i,j,k))*(vp(i,j,k)-vbr(i,j,k))
            end do
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!! For the z components of velocity.

      if(vspvar(3:3).eq.'o') then

! Damp to the GPV data.

        if(vspopt.eq.1.and.gpvvar(1:1).eq.'o') then

          do k=ksp0(1),nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              wfrc(i,j,k)=wfrc(i,j,k)-rbct(i,j,k,1)*rst8w(i,j,k)        &
     &          *(wp(i,j,k)-(wgpv(i,j,k)+wtd(i,j,k)*gtinc))
            end do
            end do

!$omp end do

          end do

! -----

! Damp to the 0.

        else

          do k=ksp0(2),nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              wfrc(i,j,k)=wfrc(i,j,k)                                   &
     &          -rbct(i,j,k,2)*rst8w(i,j,k)*wp(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

! -----

!! -----

      end if

!$omp end parallel

!!! -----

      end subroutine s_vspuvw

!-----7--------------------------------------------------------------7--

      end module m_vspuvw
