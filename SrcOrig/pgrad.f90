!***********************************************************************
      module m_pgrad
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/06/07,
!                   1999/07/05, 1999/08/03, 1999/08/23, 1999/09/30,
!                   1999/10/12, 1999/11/01, 1999/11/24, 2000/01/17,
!                   2000/12/18, 2001/02/24, 2001/06/06, 2001/06/29,
!                   2001/11/20, 2002/04/02, 2003/01/04, 2003/03/21,
!                   2003/04/30, 2003/05/19, 2003/10/31, 2003/11/05,
!                   2003/12/12, 2004/05/31, 2004/06/10, 2005/04/04,
!                   2006/06/21, 2006/11/06, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/01/30, 2009/02/27, 2011/08/09,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the pressure gradient force in the velocity equations.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_comphy
      use m_diver3d
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: pgrad, s_pgrad

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface pgrad

        module procedure s_pgrad

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
      subroutine s_pgrad(fptrnopt,fpmpopt,fpmfcopt,fpdivopt,            &
     &                   fpdx,fpdy,fpdz,fpdxiv,fpdyiv,fpdziv,           &
     &                   dts,ni,nj,nk,j31,j32,jcb,mf,mf8u,mf8v,         &
     &                   rmf,rmf8u,rmf8v,rst8u,rst8v,rst8w,u,v,wc,pp,   &
     &                   upg,vpg,wpg,tmp1,tmp2,tmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpdivopt
                       ! Formal parameter of unique index of divopt

      integer, intent(in) :: fpdx
                       ! Formal parameter of unique index of dx

      integer, intent(in) :: fpdy
                       ! Formal parameter of unique index of dy

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

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

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

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

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

! Output variables

      real, intent(out) :: upg(0:ni+1,0:nj+1,1:nk)
                       ! Pressure gradient force in u equation

      real, intent(out) :: vpg(0:ni+1,0:nj+1,1:nk)
                       ! Pressure gradient force in v equation

      real, intent(out) :: wpg(0:ni+1,0:nj+1,1:nk)
                       ! Pressure gradient force in w equation

! Internal shared variables

      integer trnopt   ! Option for terrain height setting
      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor
      integer divopt   ! Option for divergence damping

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction
      real dz          ! Grid distance in z direction

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real dziv25      ! 0.25 x dziv

      real divch       ! Divergence damping coefficient
                       ! in horizontal direction

      real divcv       ! Divergence damping coefficient
                       ! in vertical direction

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

! Remark

!     wpg: This variable is also temporary.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fptrnopt,trnopt)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpdivopt,divopt)
      call getrname(fpdx,dx)
      call getrname(fpdy,dy)
      call getrname(fpdz,dz)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      dziv25=.25e0*dziv

      divch=divndc/dts*dx*dy
      divcv=divndc/dts*dz*dz

! -----

! Calculate the divergence damping.

      if(divopt.ge.1) then

       call diver3d(idmpopt,idmfcopt,iddxiv,iddyiv,iddziv,ni,nj,nk,     &
     &              mf,rmf,rmf8u,rmf8v,rst8u,rst8v,rst8w,u,v,wc,        &
     &              tmp1,tmp2,tmp3,wpg)

      end if

! -----

!! Calculate the pressure gradient force.

!$omp parallel default(shared) private(k)

! Add the divergence damping to the pressure perturbation.

      if(divopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            tmp3(i,j,k)=pp(i,j,k)
          end do
          end do

!$omp end do

        end do

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            tmp3(i,j,k)=pp(i,j,k)+divcv*jcb(i,j,k)*tmp1(i,j,k)
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the vertical components of pressure gradient force.

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          wpg(i,j,k)=(tmp3(i,j,k-1)-tmp3(i,j,k))*dziv
        end do
        end do

!$omp end do

      end do

! -----

! Reset and add the divergence damping to the pressure perturbation for
! anisotropic case.

      if(divopt.eq.2) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            tmp3(i,j,k)=pp(i,j,k)+divch*tmp1(i,j,k)/jcb(i,j,k)
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the horizontal components of pressure gradient force.

      if(trnopt.eq.0) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            tmp3(i,j,k)=jcb(i,j,k)*tmp3(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(mfcopt.eq.0) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              upg(i,j,k)=(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vpg(i,j,k)=(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv
            end do
            end do

!$omp end do

          end do

        else

          if(mpopt.eq.0.or.mpopt.eq.10) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=mf8u(i,j)*(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv
              end do
              end do

!$omp end do

            end do

          else if(mpopt.eq.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=mf8v(i,j)*(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=mf8u(i,j)*(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=mf8v(i,j)*(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv
              end do
              end do

!$omp end do

            end do

          end if

        end if

      else if(trnopt.ge.1) then

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            tmp1(i,j,k)=((tmp3(i-1,j,k-1)+tmp3(i,j,k-1))                &
     &        +(tmp3(i-1,j,k)+tmp3(i,j,k)))*j31(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            tmp2(i,j,k)=((tmp3(i,j-1,k-1)+tmp3(i,j,k-1))                &
     &        +(tmp3(i,j-1,k)+tmp3(i,j,k)))*j32(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            tmp3(i,j,k)=jcb(i,j,k)*tmp3(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(mfcopt.eq.0) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              upg(i,j,k)=(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv               &
     &          -(tmp1(i,j,k+1)-tmp1(i,j,k))*dziv25
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vpg(i,j,k)=(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv               &
     &          -(tmp2(i,j,k+1)-tmp2(i,j,k))*dziv25
            end do
            end do

!$omp end do

          end do

        else

          if(mpopt.eq.0.or.mpopt.eq.10) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=mf8u(i,j)*((tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv  &
     &            -(tmp1(i,j,k+1)-tmp1(i,j,k))*dziv25)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=(tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv             &
     &            -(tmp2(i,j,k+1)-tmp2(i,j,k))*dziv25
              end do
              end do

!$omp end do

            end do

          else if(mpopt.eq.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=(tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv             &
     &            -(tmp1(i,j,k+1)-tmp1(i,j,k))*dziv25
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=mf8v(i,j)*((tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv  &
     &            -(tmp2(i,j,k+1)-tmp2(i,j,k))*dziv25)
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                upg(i,j,k)=mf8u(i,j)*((tmp3(i-1,j,k)-tmp3(i,j,k))*dxiv  &
     &            -(tmp1(i,j,k+1)-tmp1(i,j,k))*dziv25)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vpg(i,j,k)=mf8v(i,j)*((tmp3(i,j-1,k)-tmp3(i,j,k))*dyiv  &
     &            -(tmp2(i,j,k+1)-tmp2(i,j,k))*dziv25)
              end do
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_pgrad

!-----7--------------------------------------------------------------7--

      end module m_pgrad
