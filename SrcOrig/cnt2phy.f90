!***********************************************************************
      module m_cnt2phy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/06/01
!     Modification: 2001/06/06, 2002/04/02, 2003/01/04, 2003/03/21,
!                   2003/04/30, 2003/05/19, 2003/11/05, 2004/06/10,
!                   2006/11/06, 2007/10/19, 2008/05/02, 2008/06/09,
!                   2008/08/25, 2009/02/27, 2011/08/09, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the z components of velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: cnt2phy, s_cnt2phy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface cnt2phy

        module procedure s_cnt2phy

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
      subroutine s_cnt2phy(fpsthopt,fptrnopt,fpmpopt,fpmfcopt,          &
     &                     ni,nj,nk,j31,j32,jcb8w,mf,u,v,wc,w,          &
     &                     mf25,j31u2,j32v2)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsthopt
                       ! Formal parameter of unique index of sthopt

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

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

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

! Output variable

      real, intent(out) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

! Internal shared variables

      integer sthopt   ! Option for vertical grid stretching

      integer trnopt   ! Option for terrain height setting

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      real, intent(inout) :: mf25(0:ni+1,0:nj+1)
                       ! 0.25 x mf

      real, intent(inout) :: j31u2(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x j31 x u

      real, intent(inout) :: j32v2(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x j32 x v

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsthopt,sthopt)
      call getiname(fptrnopt,trnopt)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)

! -----

! Calculate the z components of velocity.

!$omp parallel default(shared) private(k)

      if(trnopt.eq.0) then

        if(sthopt.eq.0) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              w(i,j,k)=wc(i,j,k)
            end do
            end do

!$omp end do

          end do

        else if(sthopt.ge.1) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              w(i,j,k)=wc(i,j,k)*jcb8w(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      else

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni
            j31u2(i,j,k)=(u(i,j,k-1)+u(i,j,k))*j31(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=1,nj
          do i=1,ni-1
            j32v2(i,j,k)=(v(i,j,k-1)+v(i,j,k))*j32(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(mfcopt.eq.0) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              w(i,j,k)=wc(i,j,k)*jcb8w(i,j,k)                           &
     &          -.25e0*((j31u2(i,j,k)+j31u2(i+1,j,k))                   &
     &          +(j32v2(i,j,k)+j32v2(i,j+1,k)))
            end do
            end do

!$omp end do

          end do

        else

          if(mpopt.eq.0.or.mpopt.eq.10) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                w(i,j,k)=wc(i,j,k)*jcb8w(i,j,k)                         &
     &            -.25e0*(mf(i,j)*(j31u2(i,j,k)+j31u2(i+1,j,k))         &
     &            +(j32v2(i,j,k)+j32v2(i,j+1,k)))
              end do
              end do

!$omp end do

            end do

          else if(mpopt.eq.5) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                w(i,j,k)=wc(i,j,k)*jcb8w(i,j,k)                         &
     &            -.25e0*((j31u2(i,j,k)+j31u2(i+1,j,k))                 &
     &            +mf(i,j)*(j32v2(i,j,k)+j32v2(i,j+1,k)))
              end do
              end do

!$omp end do

            end do

          else

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              mf25(i,j)=.25e0*mf(i,j)
            end do
            end do

!$omp end do

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                w(i,j,k)=wc(i,j,k)*jcb8w(i,j,k)                         &
     &            -mf25(i,j)*((j31u2(i,j,k)+j31u2(i+1,j,k))             &
     &            +(j32v2(i,j,k)+j32v2(i,j+1,k)))
              end do
              end do

!$omp end do

            end do

          end if

        end if

      end if

!$omp end parallel

! -----

      end subroutine s_cnt2phy

!-----7--------------------------------------------------------------7--

      end module m_cnt2phy
