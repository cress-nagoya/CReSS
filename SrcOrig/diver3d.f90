!***********************************************************************
      module m_diver3d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/04/06
!     Modification: 1999/07/05, 1999/08/03, 1999/08/18, 1999/08/23,
!                   1999/10/12, 1999/11/01, 2000/01/17, 2001/02/24,
!                   2001/06/29, 2002/04/02, 2003/01/04, 2003/03/21,
!                   2003/04/30, 2003/05/19, 2004/06/10, 2006/11/06,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the 3 dimensional negative divergence.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: diver3d, s_diver3d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface diver3d

        module procedure s_diver3d

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
      subroutine s_diver3d(fpmpopt,fpmfcopt,                            &
     &                     fpdxiv,fpdyiv,fpdziv,ni,nj,nk,               &
     &                     mf,rmf,rmf8u,rmf8v,var8u,var8v,var8w,        &
     &                     u,v,wc,div3d,tmp1,tmp2,tmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

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

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: var8u(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at u points

      real, intent(in) :: var8v(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at v points

      real, intent(in) :: var8w(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at w points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

! Output variable

      real, intent(out) :: div3d(0:ni+1,0:nj+1,1:nk)
                       ! 3 dimensional divergence value

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)

! -----

!! Calculate the 3 dimensional negative divergence.

!$omp parallel default(shared) private(k)

! Optional variables at u, v and w points are multiplyed by u, v and wc.

      if(mfcopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni
            tmp1(i,j,k)=var8u(i,j,k)*u(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=1,nj
          do i=1,ni-1
            tmp2(i,j,k)=var8v(i,j,k)*v(i,j,k)
          end do
          end do

!$omp end do

        end do

      else

        if(mpopt.eq.0.or.mpopt.eq.10) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni
              tmp1(i,j,k)=var8u(i,j,k)*u(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni-1
              tmp2(i,j,k)=rmf8v(i,j,2)*var8v(i,j,k)*v(i,j,k)
            end do
            end do

!$omp end do

          end do

        else if(mpopt.eq.5) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni
              tmp1(i,j,k)=rmf8u(i,j,2)*var8u(i,j,k)*u(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni-1
              tmp2(i,j,k)=var8v(i,j,k)*v(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni
              tmp1(i,j,k)=rmf8u(i,j,2)*var8u(i,j,k)*u(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni-1
              tmp2(i,j,k)=rmf8v(i,j,2)*var8v(i,j,k)*v(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      end if

      do k=1,nk

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          tmp3(i,j,k)=var8w(i,j,k)*wc(i,j,k)
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the 3 dimensional divergence.

      if(mfcopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            div3d(i,j,k)=(tmp3(i,j,k)-tmp3(i,j,k+1))*dziv               &
     &        +((tmp1(i,j,k)-tmp1(i+1,j,k))*dxiv                        &
     &        +(tmp2(i,j,k)-tmp2(i,j+1,k))*dyiv)
          end do
          end do

!$omp end do

        end do

      else

        if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              div3d(i,j,k)=mf(i,j)*((tmp1(i,j,k)-tmp1(i+1,j,k))*dxiv    &
     &          +(tmp2(i,j,k)-tmp2(i,j+1,k))*dyiv)                      &
     &          +(tmp3(i,j,k)-tmp3(i,j,k+1))*dziv
            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              div3d(i,j,k)=rmf(i,j,1)*((tmp1(i,j,k)-tmp1(i+1,j,k))*dxiv &
     &          +(tmp2(i,j,k)-tmp2(i,j+1,k))*dyiv)                      &
     &          +(tmp3(i,j,k)-tmp3(i,j,k+1))*dziv
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_diver3d

!-----7--------------------------------------------------------------7--

      end module m_diver3d
