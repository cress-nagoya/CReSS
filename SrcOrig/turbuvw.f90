!***********************************************************************
      module m_turbuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/07/05
!     Modification: 1999/07/21, 1999/07/28, 1999/08/03, 1999/09/30,
!                   1999/10/12, 1999/11/01, 2000/01/17, 2000/12/19,
!                   2001/04/15, 2001/06/06, 2001/11/20, 2002/04/02,
!                   2003/01/04, 2003/01/20, 2003/03/13, 2003/03/21,
!                   2003/04/30, 2003/05/19, 2003/12/12, 2004/02/01,
!                   2004/06/10, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/11/06, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the velocity turbulent mixing.

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

      public :: turbuvw, s_turbuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface turbuvw

        module procedure s_turbuvw

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
      subroutine s_turbuvw(fptrnopt,fpmpopt,fpmfcopt,fpadvopt,          &
     &                     fpdxiv,fpdyiv,fpdziv,ni,nj,nk,j31,j32,       &
     &                     jcb,jcb8u,jcb8v,mf,mf8u,mf8v,rmf,rmf8u,rmf8v,&
     &                     t11,t22,t33,t12,t13,t23,t31,t32,             &
     &                     ufrc,vfrc,wfrc,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

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

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

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

! Input and output variables

      real, intent(inout) :: t11(0:ni+1,0:nj+1,1:nk)
                       ! x-x components of stress tensor

      real, intent(inout) :: t22(0:ni+1,0:nj+1,1:nk)
                       ! y-y components of stress tensor

      real, intent(inout) :: t33(0:ni+1,0:nj+1,1:nk)
                       ! z-z components of stress tensor

      real, intent(inout) :: t12(0:ni+1,0:nj+1,1:nk)
                       ! x-y components of stress tensor

      real, intent(inout) :: t13(0:ni+1,0:nj+1,1:nk)
                       ! x-z components of stress tensor

      real, intent(inout) :: t23(0:ni+1,0:nj+1,1:nk)
                       ! y-z components of stress tensor

      real, intent(inout) :: t31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of stress tensor

      real, intent(inout) :: t32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of stress tensor

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      integer trnopt   ! Option for terrain height setting
      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor
      integer advopt   ! Option for advection scheme

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real dxiv05      ! 0.5 x dxiv
      real dyiv05      ! 0.5 x dyiv
      real dxiv25      ! 0.25 x dxiv
      real dyiv25      ! 0.25 x dyiv

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

! Remark

!     t11,t22,t33,t13,t23: These variables are also temporary, because
!                          they are not used again.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fptrnopt,trnopt)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpadvopt,advopt)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      dxiv05=.5e0*dxiv
      dyiv05=.5e0*dyiv
      dxiv25=.25e0*dxiv
      dyiv25=.25e0*dyiv

! -----

!! Calculate the velocity turbulent mixing.

!$omp parallel default(shared) private(k)

! Calculate the invers of map scale factor at dot points.

      if(mfcopt.eq.1) then

        if(mpopt.eq.0.or.mpopt.eq.10) then

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-1
            t13(i,j,nk)=rmf8u(i,j-1,2)+rmf8u(i,j,2)
          end do
          end do

!$omp end do

        else if(mpopt.eq.5) then

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-1
            t23(i,j,nk)=rmf8v(i-1,j,2)+rmf8v(i,j,2)
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-1
            t13(i,j,nk)=rmf8u(i,j-1,2)+rmf8u(i,j,2)
            t23(i,j,nk)=rmf8v(i-1,j,2)+rmf8v(i,j,2)
          end do
          end do

!$omp end do

        end if

      end if

! -----

! Calculate the u turbulent mixing.

      if(trnopt.eq.0) then

        if(mfcopt.eq.0) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=1,ni-1
              t11(i,j,k)=jcb(i,j,k)*t11(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-1
              tmp1(i,j,k)=(jcb8u(i,j-1,k)+jcb8u(i,j,k))*t12(i,j,k)
            end do
            end do

!$omp end do

          end do

          if(advopt.le.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                ufrc(i,j,k)=((t11(i,j,k)-t11(i-1,j,k))*dxiv             &
     &            +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv05)                  &
     &            +(t13(i,j,k+1)-t13(i,j,k))*dziv
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                ufrc(i,j,k)=ufrc(i,j,k)+((t11(i,j,k)-t11(i-1,j,k))*dxiv &
     &            +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv05)                  &
     &            +(t13(i,j,k+1)-t13(i,j,k))*dziv
              end do
              end do

!$omp end do

            end do

          end if

        else

          if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

            if(mpopt.eq.0.or.mpopt.eq.10) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=1,ni-1
                  t11(i,j,k)=jcb(i,j,k)*t11(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-1
                  tmp1(i,j,k)=t13(i,j,nk)                               &
     &              *(jcb8u(i,j-1,k)+jcb8u(i,j,k))*t12(i,j,k)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=1,ni-1
                  t11(i,j,k)=rmf(i,j,2)*jcb(i,j,k)*t11(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-1
                  tmp1(i,j,k)=(jcb8u(i,j-1,k)+jcb8u(i,j,k))*t12(i,j,k)
                end do
                end do

!$omp end do

              end do

            end if

            if(advopt.le.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  ufrc(i,j,k)=mf8u(i,j)*((t11(i,j,k)-t11(i-1,j,k))*dxiv &
     &              +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv25)                &
     &              +(t13(i,j,k+1)-t13(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  ufrc(i,j,k)=ufrc(i,j,k)                               &
     &              +mf8u(i,j)*((t11(i,j,k)-t11(i-1,j,k))*dxiv          &
     &              +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv25)                &
     &              +(t13(i,j,k+1)-t13(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            end if

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=1,ni-1
                t11(i,j,k)=rmf(i,j,2)*jcb(i,j,k)*t11(i,j,k)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-1
                tmp1(i,j,k)=t13(i,j,nk)                                 &
     &            *(jcb8u(i,j-1,k)+jcb8u(i,j,k))*t12(i,j,k)
              end do
              end do

!$omp end do

            end do

            if(advopt.le.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  ufrc(i,j,k)                                           &
     &              =rmf8u(i,j,1)*((t11(i,j,k)-t11(i-1,j,k))*dxiv       &
     &              +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv25)                &
     &              +(t13(i,j,k+1)-t13(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  ufrc(i,j,k)=ufrc(i,j,k)                               &
     &              +rmf8u(i,j,1)*((t11(i,j,k)-t11(i-1,j,k))*dxiv       &
     &              +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv25)                &
     &              +(t13(i,j,k+1)-t13(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            end if

          end if

        end if

      else

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-1
            tmp1(i,j,k)                                                 &
     &        =(j32(i-1,j,k)+j32(i,j,k))*(t12(i,j,k-1)+t12(i,j,k))
          end do
          end do

!$omp end do

        end do

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            t13(i,j,k)=t13(i,j,k)                                       &
     &        +.125e0*(tmp1(i,j,k)+tmp1(i,j+1,k))+.25e0*j31(i,j,k)      &
     &        *((t11(i-1,j,k-1)+t11(i,j,k-1))+(t11(i-1,j,k)+t11(i,j,k)))
          end do
          end do

!$omp end do

        end do

        if(mfcopt.eq.0) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=1,ni-1
              t11(i,j,k)=jcb(i,j,k)*t11(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-1
              tmp1(i,j,k)=(jcb8u(i,j-1,k)+jcb8u(i,j,k))*t12(i,j,k)
            end do
            end do

!$omp end do

          end do

          if(advopt.le.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                ufrc(i,j,k)=((t11(i,j,k)-t11(i-1,j,k))*dxiv             &
     &            +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv05)                  &
     &            +(t13(i,j,k+1)-t13(i,j,k))*dziv
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                ufrc(i,j,k)=ufrc(i,j,k)+((t11(i,j,k)-t11(i-1,j,k))*dxiv &
     &            +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv05)                  &
     &            +(t13(i,j,k+1)-t13(i,j,k))*dziv
              end do
              end do

!$omp end do

            end do

          end if

        else

          if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

            if(mpopt.eq.0.or.mpopt.eq.10) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=1,ni-1
                  t11(i,j,k)=jcb(i,j,k)*t11(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-1
                  tmp1(i,j,k)=t13(i,j,nk)                               &
     &              *(jcb8u(i,j-1,k)+jcb8u(i,j,k))*t12(i,j,k)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=1,ni-1
                  t11(i,j,k)=rmf(i,j,2)*jcb(i,j,k)*t11(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-1
                  tmp1(i,j,k)=(jcb8u(i,j-1,k)+jcb8u(i,j,k))*t12(i,j,k)
                end do
                end do

!$omp end do

              end do

            end if

            if(advopt.le.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  ufrc(i,j,k)=mf8u(i,j)*((t11(i,j,k)-t11(i-1,j,k))*dxiv &
     &              +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv25)                &
     &              +(t13(i,j,k+1)-t13(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  ufrc(i,j,k)=ufrc(i,j,k)                               &
     &              +mf8u(i,j)*((t11(i,j,k)-t11(i-1,j,k))*dxiv          &
     &              +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv25)                &
     &              +(t13(i,j,k+1)-t13(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            end if

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=1,ni-1
                t11(i,j,k)=rmf(i,j,2)*jcb(i,j,k)*t11(i,j,k)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-1
                tmp1(i,j,k)=t13(i,j,nk)                                 &
     &            *(jcb8u(i,j-1,k)+jcb8u(i,j,k))*t12(i,j,k)
              end do
              end do

!$omp end do

            end do

            if(advopt.le.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  ufrc(i,j,k)                                           &
     &              =rmf8u(i,j,1)*((t11(i,j,k)-t11(i-1,j,k))*dxiv       &
     &              +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv25)                &
     &              +(t13(i,j,k+1)-t13(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  ufrc(i,j,k)=ufrc(i,j,k)                               &
     &              +rmf8u(i,j,1)*((t11(i,j,k)-t11(i-1,j,k))*dxiv       &
     &              +(tmp1(i,j+1,k)-tmp1(i,j,k))*dyiv25)                &
     &              +(t13(i,j,k+1)-t13(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            end if

          end if

        end if

      end if

! -----

! Calculate the v turbulent mixing.

      if(trnopt.eq.0) then

        if(mfcopt.eq.0) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=2,ni-2
              t22(i,j,k)=jcb(i,j,k)*t22(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-1
              tmp1(i,j,k)=(jcb8v(i-1,j,k)+jcb8v(i,j,k))*t12(i,j,k)
            end do
            end do

!$omp end do

          end do

          if(advopt.le.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vfrc(i,j,k)=((t22(i,j,k)-t22(i,j-1,k))*dyiv             &
     &            +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv05)                  &
     &            +(t23(i,j,k+1)-t23(i,j,k))*dziv
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vfrc(i,j,k)=vfrc(i,j,k)+((t22(i,j,k)-t22(i,j-1,k))*dyiv &
     &            +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv05)                  &
     &            +(t23(i,j,k+1)-t23(i,j,k))*dziv
              end do
              end do

!$omp end do

            end do

          end if

        else

          if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

            if(mpopt.eq.0.or.mpopt.eq.10) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=2,ni-2
                  t22(i,j,k)=rmf(i,j,2)*jcb(i,j,k)*t22(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-1
                  tmp1(i,j,k)=(jcb8v(i-1,j,k)+jcb8v(i,j,k))*t12(i,j,k)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=2,ni-2
                  t22(i,j,k)=jcb(i,j,k)*t22(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-1
                  tmp1(i,j,k)=t23(i,j,nk)                               &
     &              *(jcb8v(i-1,j,k)+jcb8v(i,j,k))*t12(i,j,k)
                end do
                end do

!$omp end do

              end do

            end if

            if(advopt.le.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  vfrc(i,j,k)=mf8v(i,j)*((t22(i,j,k)-t22(i,j-1,k))*dyiv &
     &              +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv25)                &
     &              +(t23(i,j,k+1)-t23(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  vfrc(i,j,k)=vfrc(i,j,k)                               &
     &              +mf8v(i,j)*((t22(i,j,k)-t22(i,j-1,k))*dyiv          &
     &              +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv25)                &
     &              +(t23(i,j,k+1)-t23(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            end if

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=2,ni-2
                t22(i,j,k)=rmf(i,j,2)*jcb(i,j,k)*t22(i,j,k)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-1
                tmp1(i,j,k)=t23(i,j,nk)                                 &
     &            *(jcb8v(i-1,j,k)+jcb8v(i,j,k))*t12(i,j,k)
              end do
              end do

!$omp end do

            end do

            if(advopt.le.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  vfrc(i,j,k)                                           &
     &              =rmf8v(i,j,1)*((t22(i,j,k)-t22(i,j-1,k))*dyiv       &
     &              +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv25)                &
     &              +(t23(i,j,k+1)-t23(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  vfrc(i,j,k)=vfrc(i,j,k)                               &
     &              +rmf8v(i,j,1)*((t22(i,j,k)-t22(i,j-1,k))*dyiv       &
     &              +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv25)                &
     &              +(t23(i,j,k+1)-t23(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            end if

          end if

        end if

      else

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-1
            tmp1(i,j,k)                                                 &
     &        =(j31(i,j-1,k)+j31(i,j,k))*(t12(i,j,k-1)+t12(i,j,k))
          end do
          end do

!$omp end do

        end do

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            t23(i,j,k)=t23(i,j,k)                                       &
     &        +.125e0*(tmp1(i,j,k)+tmp1(i+1,j,k))+.25e0*j32(i,j,k)      &
     &        *((t22(i,j-1,k-1)+t22(i,j,k-1))+(t22(i,j-1,k)+t22(i,j,k)))
          end do
          end do

!$omp end do

        end do

        if(mfcopt.eq.0) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=2,ni-2
              t22(i,j,k)=jcb(i,j,k)*t22(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-1
              tmp1(i,j,k)=(jcb8v(i-1,j,k)+jcb8v(i,j,k))*t12(i,j,k)
            end do
            end do

!$omp end do

          end do

          if(advopt.le.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vfrc(i,j,k)=((t22(i,j,k)-t22(i,j-1,k))*dyiv             &
     &            +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv05)                  &
     &            +(t23(i,j,k+1)-t23(i,j,k))*dziv
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                vfrc(i,j,k)=vfrc(i,j,k)+((t22(i,j,k)-t22(i,j-1,k))*dyiv &
     &            +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv05)                  &
     &            +(t23(i,j,k+1)-t23(i,j,k))*dziv
              end do
              end do

!$omp end do

            end do

          end if

        else

          if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

            if(mpopt.eq.0.or.mpopt.eq.10) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=2,ni-2
                  t22(i,j,k)=rmf(i,j,2)*jcb(i,j,k)*t22(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-1
                  tmp1(i,j,k)=(jcb8v(i-1,j,k)+jcb8v(i,j,k))*t12(i,j,k)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=2,ni-2
                  t22(i,j,k)=jcb(i,j,k)*t22(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-1
                  tmp1(i,j,k)=t23(i,j,nk)                               &
     &              *(jcb8v(i-1,j,k)+jcb8v(i,j,k))*t12(i,j,k)
                end do
                end do

!$omp end do

              end do

            end if

            if(advopt.le.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  vfrc(i,j,k)=mf8v(i,j)*((t22(i,j,k)-t22(i,j-1,k))*dyiv &
     &              +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv25)                &
     &              +(t23(i,j,k+1)-t23(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  vfrc(i,j,k)=vfrc(i,j,k)                               &
     &              +mf8v(i,j)*((t22(i,j,k)-t22(i,j-1,k))*dyiv          &
     &              +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv25)                &
     &              +(t23(i,j,k+1)-t23(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            end if

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=2,ni-2
                t22(i,j,k)=rmf(i,j,2)*jcb(i,j,k)*t22(i,j,k)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-1
                tmp1(i,j,k)=t23(i,j,nk)                                 &
     &            *(jcb8v(i-1,j,k)+jcb8v(i,j,k))*t12(i,j,k)
              end do
              end do

!$omp end do

            end do

            if(advopt.le.3) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  vfrc(i,j,k)                                           &
     &              =rmf8v(i,j,1)*((t22(i,j,k)-t22(i,j-1,k))*dyiv       &
     &              +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv25)                &
     &              +(t23(i,j,k+1)-t23(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  vfrc(i,j,k)=vfrc(i,j,k)                               &
     &              +rmf8v(i,j,1)*((t22(i,j,k)-t22(i,j-1,k))*dyiv       &
     &              +(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv25)                &
     &              +(t23(i,j,k+1)-t23(i,j,k))*dziv
                end do
                end do

!$omp end do

              end do

            end if

          end if

        end if

      end if

! -----

! Calculate the w turbulent mixing.

      if(trnopt.eq.0) then

        if(mfcopt.eq.0) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              t11(i,j,k)=(jcb8u(i,j,k-1)+jcb8u(i,j,k))*t31(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              t22(i,j,k)=(jcb8v(i,j,k-1)+jcb8v(i,j,k))*t32(i,j,k)
            end do
            end do

!$omp end do

          end do

          if(advopt.le.3) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                wfrc(i,j,k)=(t33(i,j,k)-t33(i,j,k-1))*dziv              &
     &            +((t11(i+1,j,k)-t11(i,j,k))*dxiv05                    &
     &            +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                wfrc(i,j,k)=wfrc(i,j,k)+(t33(i,j,k)-t33(i,j,k-1))*dziv  &
     &            +((t11(i+1,j,k)-t11(i,j,k))*dxiv05                    &
     &            +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
              end do
              end do

!$omp end do

            end do

          end if

        else

          if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

            if(mpopt.eq.0.or.mpopt.eq.10) then

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  t11(i,j,k)=(jcb8u(i,j,k-1)+jcb8u(i,j,k))*t31(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  t22(i,j,k)=rmf8v(i,j,2)                               &
     &              *(jcb8v(i,j,k-1)+jcb8v(i,j,k))*t32(i,j,k)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  t11(i,j,k)=rmf8u(i,j,2)                               &
     &              *(jcb8u(i,j,k-1)+jcb8u(i,j,k))*t31(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  t22(i,j,k)=(jcb8v(i,j,k-1)+jcb8v(i,j,k))*t32(i,j,k)
                end do
                end do

!$omp end do

              end do

            end if

            if(advopt.le.3) then

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  wfrc(i,j,k)=(t33(i,j,k)-t33(i,j,k-1))*dziv            &
     &              +mf(i,j)*((t11(i+1,j,k)-t11(i,j,k))*dxiv05          &
     &              +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  wfrc(i,j,k)=wfrc(i,j,k)+(t33(i,j,k)-t33(i,j,k-1))*dziv&
     &              +mf(i,j)*((t11(i+1,j,k)-t11(i,j,k))*dxiv05          &
     &              +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
                end do
                end do

!$omp end do

              end do

            end if

          else

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                t11(i,j,k)=rmf8u(i,j,2)                                 &
     &            *(jcb8u(i,j,k-1)+jcb8u(i,j,k))*t31(i,j,k)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                t22(i,j,k)=rmf8v(i,j,2)                                 &
     &            *(jcb8v(i,j,k-1)+jcb8v(i,j,k))*t32(i,j,k)
              end do
              end do

!$omp end do

            end do

            if(advopt.le.3) then

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  wfrc(i,j,k)=(t33(i,j,k)-t33(i,j,k-1))*dziv            &
     &              +rmf(i,j,1)*((t11(i+1,j,k)-t11(i,j,k))*dxiv05       &
     &              +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  wfrc(i,j,k)=wfrc(i,j,k)+(t33(i,j,k)-t33(i,j,k-1))*dziv&
     &              +rmf(i,j,1)*((t11(i+1,j,k)-t11(i,j,k))*dxiv05       &
     &              +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
                end do
                end do

!$omp end do

              end do

            end if

          end if

        end if

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            t11(i,j,k)                                                  &
     &        =(j31(i,j,k)+j31(i,j,k+1))*(t31(i,j,k)+t31(i,j,k+1))
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            t22(i,j,k)                                                  &
     &        =(j32(i,j,k)+j32(i,j,k+1))*(t32(i,j,k)+t32(i,j,k+1))
          end do
          end do

!$omp end do

        end do

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            t33(i,j,k)=t33(i,j,k)+.125e0*((t11(i,j,k)+t11(i+1,j,k))     &
     &        +(t22(i,j,k)+t22(i,j+1,k)))
          end do
          end do

!$omp end do

        end do

        if(mfcopt.eq.0) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              t11(i,j,k)=(jcb8u(i,j,k-1)+jcb8u(i,j,k))*t31(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              t22(i,j,k)=(jcb8v(i,j,k-1)+jcb8v(i,j,k))*t32(i,j,k)
            end do
            end do

!$omp end do

          end do

          if(advopt.le.3) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                wfrc(i,j,k)=(t33(i,j,k)-t33(i,j,k-1))*dziv              &
     &            +((t11(i+1,j,k)-t11(i,j,k))*dxiv05                    &
     &            +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                wfrc(i,j,k)=wfrc(i,j,k)+(t33(i,j,k)-t33(i,j,k-1))*dziv  &
     &            +((t11(i+1,j,k)-t11(i,j,k))*dxiv05                    &
     &            +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
              end do
              end do

!$omp end do

            end do

          end if

        else

          if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

            if(mpopt.eq.0.or.mpopt.eq.10) then

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  t11(i,j,k)=(jcb8u(i,j,k-1)+jcb8u(i,j,k))*t31(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  t22(i,j,k)=rmf8v(i,j,2)                               &
     &              *(jcb8v(i,j,k-1)+jcb8v(i,j,k))*t32(i,j,k)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  t11(i,j,k)=rmf8u(i,j,2)                               &
     &              *(jcb8u(i,j,k-1)+jcb8u(i,j,k))*t31(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  t22(i,j,k)=(jcb8v(i,j,k-1)+jcb8v(i,j,k))*t32(i,j,k)
                end do
                end do

!$omp end do

              end do

            end if

            if(advopt.le.3) then

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  wfrc(i,j,k)=(t33(i,j,k)-t33(i,j,k-1))*dziv            &
     &              +mf(i,j)*((t11(i+1,j,k)-t11(i,j,k))*dxiv05          &
     &              +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  wfrc(i,j,k)=wfrc(i,j,k)+(t33(i,j,k)-t33(i,j,k-1))*dziv&
     &              +mf(i,j)*((t11(i+1,j,k)-t11(i,j,k))*dxiv05          &
     &              +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
                end do
                end do

!$omp end do

              end do

            end if

          else

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                t11(i,j,k)=rmf8u(i,j,2)                                 &
     &            *(jcb8u(i,j,k-1)+jcb8u(i,j,k))*t31(i,j,k)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                t22(i,j,k)=rmf8v(i,j,2)                                 &
     &            *(jcb8v(i,j,k-1)+jcb8v(i,j,k))*t32(i,j,k)
              end do
              end do

!$omp end do

            end do

            if(advopt.le.3) then

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  wfrc(i,j,k)=(t33(i,j,k)-t33(i,j,k-1))*dziv            &
     &              +rmf(i,j,1)*((t11(i+1,j,k)-t11(i,j,k))*dxiv05       &
     &              +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-2
                  wfrc(i,j,k)=wfrc(i,j,k)+(t33(i,j,k)-t33(i,j,k-1))*dziv&
     &              +rmf(i,j,1)*((t11(i+1,j,k)-t11(i,j,k))*dxiv05       &
     &              +(t22(i,j+1,k)-t22(i,j,k))*dyiv05)
                end do
                end do

!$omp end do

              end do

            end if

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_turbuvw

!-----7--------------------------------------------------------------7--

      end module m_turbuvw
