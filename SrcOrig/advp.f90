!***********************************************************************
      module m_advp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/08/18, 1999/09/30, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/04/18, 2000/12/18, 2001/04/15,
!                   2001/05/29, 2001/06/06, 2001/11/20, 2002/04/02,
!                   2003/01/04, 2003/03/21, 2003/04/30, 2003/05/19,
!                   2003/11/28, 2004/03/05, 2004/04/15, 2004/06/10,
!                   2004/08/01, 2006/02/13, 2006/04/03, 2006/05/12,
!                   2006/11/06, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the pressure advection.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: advp, s_advp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface advp

        module procedure s_advp

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
      subroutine s_advp(fpadvopt,fpmpopt,fpmfcopt,                      &
     &                  fpdiaopt,fpiwest,fpieast,fpjsouth,fpjnorth,     &
     &                  fpdxiv,fpdyiv,fpdziv,ni,nj,nk,mf8u,mf8v,        &
     &                  jcb8u,jcb8v,jcb8w,u,v,wc,pp,pfrc,               &
     &                  jcbxu,jcbxv,jcbxwc,hadv,vadv,tmp1,tmp2,tmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpdiaopt
                       ! Formal parameter of unique index of diaopt

      integer, intent(in) :: fpiwest
                       ! Formal parameter of unique index of iwest

      integer, intent(in) :: fpieast
                       ! Formal parameter of unique index of ieast

      integer, intent(in) :: fpjsouth
                       ! Formal parameter of unique index of jsouth

      integer, intent(in) :: fpjnorth
                       ! Formal parameter of unique index of jnorth

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

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

! Input and output variable

      real, intent(inout) :: pfrc(0:ni+1,0:nj+1,1:nk)
                       ! Pressure forcing term

! Internal shared variables

      integer advopt   ! Option for advection scheme
      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor
      integer diaopt   ! Option for diabatic calculation

      integer iwest    ! Added index on west boundary
      integer jsouth   ! Added index on south boundary

      integer ieast    ! Subtracted index on east boundary
      integer jnorth   ! Subtracted index on north boundary

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real dxv05n      ! - 0.5 x dxiv
      real dyv05n      ! - 0.5 x dyiv
      real dzv05n      ! - 0.5 x dziv

      real dxv24       ! dxiv / 24.0
      real dyv24       ! dyiv / 24.0
      real dzv24       ! dziv / 24.0

      real, intent(inout) :: jcbxu(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points x u

      real, intent(inout) :: jcbxv(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points x v

      real, intent(inout) :: jcbxwc(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points x wc

      real, intent(inout) :: hadv(0:ni+1,0:nj+1,1:nk)
                       ! Advection value in x and y direction

      real, intent(inout) :: vadv(0:ni+1,0:nj+1,1:nk)
                       ! Advection value in z direction

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

      call getiname(fpadvopt,advopt)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpdiaopt,diaopt)
      call getiname(fpiwest,iwest)
      call getiname(fpieast,ieast)
      call getiname(fpjsouth,jsouth)
      call getiname(fpjnorth,jnorth)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      dxv05n=-.5e0*dxiv
      dyv05n=-.5e0*dyiv
      dzv05n=-.5e0*dziv

      dxv24=oned24*dxiv
      dyv24=oned24*dyiv
      dzv24=oned24*dziv

! -----

!!! Calculate the pressure advection.

!$omp parallel default(shared) private(k)

!! Perform the centered fdm scheme.

      if(advopt.le.3) then

! The variables at the u, v and w points is multiplyed by u, v and w.

        if(mfcopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni
              jcbxu(i,j,k)=jcb8u(i,j,k)*u(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1,nj
            do i=1,ni-1
              jcbxv(i,j,k)=jcb8v(i,j,k)*v(i,j,k)
            end do
            end do

!$omp end do

          end do

          do k=1,nk

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              jcbxwc(i,j,k)=jcb8w(i,j,k)*wc(i,j,k)
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
                jcbxu(i,j,k)=mf8u(i,j)*jcb8u(i,j,k)*u(i,j,k)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=1,nj
              do i=1,ni-1
                jcbxv(i,j,k)=jcb8v(i,j,k)*v(i,j,k)
              end do
              end do

!$omp end do

            end do

          else if(mpopt.eq.5) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni
                jcbxu(i,j,k)=jcb8u(i,j,k)*u(i,j,k)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=1,nj
              do i=1,ni-1
                jcbxv(i,j,k)=mf8v(i,j)*jcb8v(i,j,k)*v(i,j,k)
              end do
              end do

!$omp end do

            end do

          else

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni
                jcbxu(i,j,k)=mf8u(i,j)*jcb8u(i,j,k)*u(i,j,k)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=1,nj
              do i=1,ni-1
                jcbxv(i,j,k)=mf8v(i,j)*jcb8v(i,j,k)*v(i,j,k)
              end do
              end do

!$omp end do

            end do

          end if

          do k=1,nk

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              jcbxwc(i,j,k)=jcb8w(i,j,k)*wc(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

! -----

! Calculate the 2nd order pressure advection.

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            tmp1(i,j,k)=jcbxu(i,j,k)*(pp(i,j,k)-pp(i-1,j,k))*dxv05n
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            tmp2(i,j,k)=jcbxv(i,j,k)*(pp(i,j,k)-pp(i,j-1,k))*dyv05n
          end do
          end do

!$omp end do

        end do

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            tmp3(i,j,k)=jcbxwc(i,j,k)*(pp(i,j,k)-pp(i,j,k-1))*dzv05n
          end do
          end do

!$omp end do

        end do

        if(diaopt.eq.0) then

          if(advopt.eq.1) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                pfrc(i,j,k)=(tmp3(i,j,k)+tmp3(i,j,k+1))                 &
     &            +((tmp1(i,j,k)+tmp1(i+1,j,k))                         &
     &            +(tmp2(i,j,k)+tmp2(i,j+1,k)))
              end do
              end do

!$omp end do

            end do

          else

            if(advopt.eq.2) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  pfrc(i,j,k)=(tmp1(i,j,k)+tmp1(i+1,j,k))               &
     &              +(tmp2(i,j,k)+tmp2(i,j+1,k))
                end do
                end do

!$omp end do

              end do

            else if(advopt.eq.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  pfrc(i,j,k)=(tmp1(i,j,k)+tmp1(i+1,j,k))               &
     &              +(tmp2(i,j,k)+tmp2(i,j+1,k))

                  vadv(i,j,k)=tmp3(i,j,k)+tmp3(i,j,k+1)

                end do
                end do

!$omp end do

              end do

            end if

          end if

        else

          if(advopt.eq.1) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                pfrc(i,j,k)=pfrc(i,j,k)+((tmp3(i,j,k)+tmp3(i,j,k+1))    &
     &            +((tmp1(i,j,k)+tmp1(i+1,j,k))                         &
     &            +(tmp2(i,j,k)+tmp2(i,j+1,k))))
              end do
              end do

!$omp end do

            end do

          else

            if(advopt.eq.2) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  hadv(i,j,k)=(tmp1(i,j,k)+tmp1(i+1,j,k))               &
     &              +(tmp2(i,j,k)+tmp2(i,j+1,k))
                end do
                end do

!$omp end do

              end do

            else if(advopt.eq.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  hadv(i,j,k)=(tmp1(i,j,k)+tmp1(i+1,j,k))               &
     &              +(tmp2(i,j,k)+tmp2(i,j+1,k))

                  vadv(i,j,k)=tmp3(i,j,k)+tmp3(i,j,k+1)

                end do
                end do

!$omp end do

              end do

            end if

          end if

        end if

! -----

! Calculate the 4th order pressure advection.

        if(advopt.eq.2.or.advopt.eq.3) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2+jsouth,nj-2-jnorth
            do i=1+iwest,ni-1-ieast
              tmp1(i,j,k)=(jcbxu(i,j,k)+jcbxu(i+1,j,k))                 &
     &          *(pp(i+1,j,k)-pp(i-1,j,k))*dxv24
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-1-jnorth
            do i=2+iwest,ni-2-ieast
              tmp2(i,j,k)=(jcbxv(i,j,k)+jcbxv(i,j+1,k))                 &
     &          *(pp(i,j+1,k)-pp(i,j-1,k))*dyv24
            end do
            end do

!$omp end do

          end do

          if(diaopt.eq.0) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2+jsouth,nj-2-jnorth
              do i=2+iwest,ni-2-ieast
                pfrc(i,j,k)=fourd3*pfrc(i,j,k)                          &
     &            +((tmp1(i-1,j,k)+tmp1(i+1,j,k))                       &
     &            +(tmp2(i,j-1,k)+tmp2(i,j+1,k)))
              end do
              end do

!$omp end do

            end do

            if(advopt.eq.2) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  pfrc(i,j,k)=pfrc(i,j,k)+(tmp3(i,j,k)+tmp3(i,j,k+1))
                end do
                end do

!$omp end do

              end do

            else if(advopt.eq.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2+jsouth,nj-2-jnorth
                do i=2+iwest,ni-2-ieast
                  tmp3(i,j,k)=(jcbxwc(i,j,k)+jcbxwc(i,j,k+1))           &
     &              *(pp(i,j,k+1)-pp(i,j,k-1))*dzv24
                end do
                end do

!$omp end do

              end do

              do k=3,nk-3

!$omp do schedule(runtime) private(i,j)

                do j=2+jsouth,nj-2-jnorth
                do i=2+iwest,ni-2-ieast
                  vadv(i,j,k)                                           &
     &              =fourd3*vadv(i,j,k)+(tmp3(i,j,k-1)+tmp3(i,j,k+1))
                end do
                end do

!$omp end do

              end do

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  pfrc(i,j,k)=pfrc(i,j,k)+vadv(i,j,k)
                end do
                end do

!$omp end do

              end do

            end if

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2+jsouth,nj-2-jnorth
              do i=2+iwest,ni-2-ieast
                hadv(i,j,k)=fourd3*hadv(i,j,k)                          &
     &            +((tmp1(i-1,j,k)+tmp1(i+1,j,k))                       &
     &            +(tmp2(i,j-1,k)+tmp2(i,j+1,k)))
              end do
              end do

!$omp end do

            end do

            if(advopt.eq.2) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  pfrc(i,j,k)=pfrc(i,j,k)                               &
     &              +(hadv(i,j,k)+(tmp3(i,j,k)+tmp3(i,j,k+1)))
                end do
                end do

!$omp end do

              end do

            else if(advopt.eq.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2+jsouth,nj-2-jnorth
                do i=2+iwest,ni-2-ieast
                  tmp3(i,j,k)=(jcbxwc(i,j,k)+jcbxwc(i,j,k+1))           &
     &              *(pp(i,j,k+1)-pp(i,j,k-1))*dzv24
                end do
                end do

!$omp end do

              end do

              do k=3,nk-3

!$omp do schedule(runtime) private(i,j)

                do j=2+jsouth,nj-2-jnorth
                do i=2+iwest,ni-2-ieast
                  vadv(i,j,k)                                           &
     &              =fourd3*vadv(i,j,k)+(tmp3(i,j,k-1)+tmp3(i,j,k+1))
                end do
                end do

!$omp end do

              end do

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  pfrc(i,j,k)=pfrc(i,j,k)+(hadv(i,j,k)+vadv(i,j,k))
                end do
                end do

!$omp end do

              end do

            end if

          end if

        end if

! -----

!! -----

!! Set the forcing term with 0 in the case the Cubic Lagrange scheme is
!! performed.

      else

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            pfrc(i,j,k)=0.e0
          end do
          end do

!$omp end do

        end do

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_advp

!-----7--------------------------------------------------------------7--

      end module m_advp
