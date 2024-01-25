!***********************************************************************
      module m_rstuvwc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/01/04
!     Modification: 2003/03/21, 2003/04/30, 2003/05/19, 2004/06/10,
!                   2004/09/10, 2006/11/06, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/08/09, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     the base state density x the Jacobian is multiplyed by u, v and w.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rstuvwc, s_rstuvwc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rstuvwc

        module procedure s_rstuvwc

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
      subroutine s_rstuvwc(fpmpopt,fpmfcopt,                            &
     &                     fpiwest,fpieast,fpjsouth,fpjnorth,           &
     &                     ni,nj,nk,mf8u,mf8v,rst8u,rst8v,rst8w,        &
     &                     u,v,wc,rstxu,rstxv,rstxwc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpiwest
                       ! Formal parameter of unique index of iwest

      integer, intent(in) :: fpieast
                       ! Formal parameter of unique index of ieast

      integer, intent(in) :: fpjsouth
                       ! Formal parameter of unique index of jsouth

      integer, intent(in) :: fpjnorth
                       ! Formal parameter of unique index of jnorth

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(in) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

! Output variables

      real, intent(out) :: rstxu(0:ni+1,0:nj+1,1:nk)
                       ! u x base state density x Jacobian

      real, intent(out) :: rstxv(0:ni+1,0:nj+1,1:nk)
                       ! v x base state density x Jacobian

      real, intent(out) :: rstxwc(0:ni+1,0:nj+1,1:nk)
                       ! wc x base state density x Jacobian

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      integer iwest    ! Added index on west boundary
      integer jsouth   ! Added index on south boundary

      integer ieast    ! Subtracted index on east boundary
      integer jnorth   ! Subtracted index on north boundary

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpiwest,iwest)
      call getiname(fpieast,ieast)
      call getiname(fpjsouth,jsouth)
      call getiname(fpjnorth,jnorth)

! -----

! The base state density x the Jacobian is multiplyed by u, v and w.

!$omp parallel default(shared) private(k)

      if(mfcopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=iwest,ni+1-ieast
            rstxu(i,j,k)=rst8u(i,j,k)*u(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=jsouth,nj+1-jnorth
          do i=1,ni-1
            rstxv(i,j,k)=rst8v(i,j,k)*v(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            rstxwc(i,j,k)=rst8w(i,j,k)*wc(i,j,k)
          end do
          end do

!$omp end do

        end do

      else

        if(mpopt.eq.0.or.mpopt.eq.10) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=iwest,ni+1-ieast
              rstxu(i,j,k)=mf8u(i,j)*rst8u(i,j,k)*u(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=jsouth,nj+1-jnorth
            do i=1,ni-1
              rstxv(i,j,k)=rst8v(i,j,k)*v(i,j,k)
            end do
            end do

!$omp end do

          end do

        else if(mpopt.eq.5) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=iwest,ni+1-ieast
              rstxu(i,j,k)=rst8u(i,j,k)*u(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=jsouth,nj+1-jnorth
            do i=1,ni-1
              rstxv(i,j,k)=mf8v(i,j)*rst8v(i,j,k)*v(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=iwest,ni+1-ieast
              rstxu(i,j,k)=mf8u(i,j)*rst8u(i,j,k)*u(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=jsouth,nj+1-jnorth
            do i=1,ni-1
              rstxv(i,j,k)=mf8v(i,j)*rst8v(i,j,k)*v(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            rstxwc(i,j,k)=rst8w(i,j,k)*wc(i,j,k)
          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! -----

      end subroutine s_rstuvwc

!-----7--------------------------------------------------------------7--

      end module m_rstuvwc
