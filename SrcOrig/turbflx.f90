!***********************************************************************
      module m_turbflx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/07/05
!     Modification: 1999/07/21, 1999/08/18, 1999/09/30, 1999/10/12,
!                   1999/11/01, 2000/01/17, 2000/12/19, 2001/02/24,
!                   2001/06/06, 2001/11/20, 2002/01/15, 2002/04/02,
!                   2003/01/04, 2003/03/21, 2003/04/30, 2003/05/19,
!                   2003/12/12, 2004/02/01, 2004/06/10, 2006/02/13,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the turbulent fluxes for optional scalar variable.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc8w
      use m_comindx
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: turbflx, s_turbflx

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface turbflx

        module procedure s_turbflx

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
      subroutine s_turbflx(fptrnopt,fpsfcopt,fpdxiv,fpdyiv,fpdziv,      &
     &                     ni,nj,nk,j31,j32,jcb,s,sfrc,rkh8u,rkh8v,     &
     &                     rkv8w,h1,h2,h3,j31s,j32s,jcbs)
!***********************************************************************

! Input variables

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

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

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: s(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable

      real, intent(in) :: sfrc(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar forcing term

      real, intent(in) :: rkh8u(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity
                       ! / Jacobian at u points

      real, intent(in) :: rkh8v(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity
                       ! / Jacobian at v points

      real, intent(in) :: rkv8w(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy diffusivity
                       ! / Jacobian at w points

! Output variables

      real, intent(out) :: h1(0:ni+1,0:nj+1,1:nk)
                       ! x components of turbulent fluxes

      real, intent(out) :: h2(0:ni+1,0:nj+1,1:nk)
                       ! y components of turbulent fluxes

      real, intent(out) :: h3(0:ni+1,0:nj+1,1:nk)
                       ! z components of turbulent fluxes

! Internal shared variables

      integer trnopt   ! Option for terrain height setting
      integer sfcopt   ! Option for surface physics

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real dxv05       ! 0.5 x dxiv
      real dyv05       ! 0.5 x dyiv

      real dzv125      ! 0.125 x dziv

      real, intent(inout) :: j31s(0:ni+1,0:nj+1,1:nk)
                       ! 4.0 x j31 x optional scalar variable

      real, intent(inout) :: j32s(0:ni+1,0:nj+1,1:nk)
                       ! 4.0 x j31 x optional scalar variable

      real, intent(inout) :: jcbs(0:ni+1,0:nj+1,1:nk)
                       ! jcb x optional scalar variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the namelist variables.

      call getiname(fptrnopt,trnopt)
      call getiname(fpsfcopt,sfcopt)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      dxv05=.5e0*dxiv
      dyv05=.5e0*dyiv

      dzv125=.125e0*dziv

! -----

! Set the common used array.

      if(trnopt.ge.1) then

!$omp parallel default(shared) private(k)

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-1
            j31s(i,j,k)=((s(i-1,j,k-1)+s(i,j,k-1))                      &
     &        +(s(i-1,j,k)+s(i,j,k)))*j31(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=1,ni-1
            j32s(i,j,k)=((s(i,j-1,k-1)+s(i,j,k-1))                      &
     &        +(s(i,j-1,k)+s(i,j,k)))*j32(i,j,k)
          end do
          end do

!$omp end do

        end do

!$omp end parallel

        call bc8w(idbbc,idtbc,ni,nj,nk,j31s)
        call bc8w(idbbc,idtbc,ni,nj,nk,j32s)

      end if

! -----

!! Calculate the turbulent fluxes.

!$omp parallel default(shared) private(k)

! Set the common used array.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          jcbs(i,j,k)=s(i,j,k)*jcb(i,j,k)
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the x components of the turbulent fluxes.

      if(trnopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            h1(i,j,k)=rkh8u(i,j,k)*(jcbs(i,j,k)-jcbs(i-1,j,k))*dxv05
          end do
          end do

!$omp end do

        end do

      else if(trnopt.ge.1) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            h1(i,j,k)=rkh8u(i,j,k)*((jcbs(i,j,k)-jcbs(i-1,j,k))*dxv05   &
     &        +(j31s(i,j,k+1)-j31s(i,j,k))*dzv125)
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the y components of the turbulent fluxes.

      if(trnopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            h2(i,j,k)=rkh8v(i,j,k)*(jcbs(i,j,k)-jcbs(i,j-1,k))*dyv05
          end do
          end do

!$omp end do

        end do

      else if(trnopt.ge.1) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            h2(i,j,k)=rkh8v(i,j,k)*((jcbs(i,j,k)-jcbs(i,j-1,k))*dyv05   &
     &        +(j32s(i,j,k+1)-j32s(i,j,k))*dzv125)
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Get the z components of the turbulent fluxes.

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          h3(i,j,k)=rkv8w(i,j,k)*(s(i,j,k)-s(i,j,k-1))*dziv
        end do
        end do

!$omp end do

      end do

      if(sfcopt.ge.1) then

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          h3(i,j,2)=sfrc(i,j,1)
        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_turbflx

!-----7--------------------------------------------------------------7--

      end module m_turbflx
