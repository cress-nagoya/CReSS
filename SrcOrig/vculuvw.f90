!***********************************************************************
      module m_vculuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/04/03
!     Modification: 2006/05/12, 2006/06/21, 2006/09/30, 2007/01/05,
!                   2007/07/30, 2007/10/19, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the horizontal velocity advection by the Cubic Lagrange
!     scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_copy3d
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vculuvw, s_vculuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vculuvw

        module procedure s_vculuvw

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
      subroutine s_vculuvw(fpdziv,ivstp,dtsep,ni,nj,nk,wc,wc8s,up,vp,wp,&
     &                     uf,vf,wf)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: ivstp
                       ! Index of time steps of
                       ! vertical Cubic Lagrange advection

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtsep
                       ! Time steps interval of
                       ! vertical Cubic Lagrange advection

      real, intent(in) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

      real, intent(in) :: wc8s(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity
                       ! at scalar points

! Input and output variables

      real, intent(inout) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(inout) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(inout) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(inout) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(inout) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(inout) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

! Internal shared variables

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2
      integer nkm3     ! nk - 3
      integer nkm4     ! nk - 4

      real dziv        ! Inverse of dz

      real dzt22       ! dziv^2 x dtsep^2 / 2.0
      real dzt36       ! - dziv^3 x dtsep^3 / 6.0
      real dzt16       ! - dziv x dtsep / 6.0
      real dzt1        ! - dziv x dtsep

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real wc8u        ! wc at u points
      real wc8v        ! wc at v points

      real a           ! Temporary variable
      real b           ! Temporary variable
      real c           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2
      nkm3=nk-3
      nkm4=nk-4

      dzt22=.5e0*dziv*dziv*dtsep*dtsep
      dzt36=-oned6*dziv*dziv*dziv*dtsep*dtsep*dtsep
      dzt16=-oned6*dziv*dtsep
      dzt1=-dziv*dtsep

! -----

! Copy the future variable to the past.

      if(ivstp.ge.2) then

        call copy3d(0,ni+1,0,nj+1,1,nk,uf,up)
        call copy3d(0,ni+1,0,nj+1,1,nk,vf,vp)
        call copy3d(0,ni+1,0,nj+1,1,nk,wf,wp)

      end if

! -----

!! Calculate the velocity advection vertically.

!$omp parallel default(shared) private(k)

! Calculate the u advection vertically.

      do k=3,nk-3

!$omp do schedule(runtime) private(i,j,wc8u,a,b,c)

        do j=2,nj-2
        do i=2,ni-1

          wc8u=.5e0*(wc8s(i-1,j,k)+wc8s(i,j,k))

          if(wc8u.gt.0.e0) then

            a=dzt36*(up(i,j,k+1)                                        &
     &        -3.e0*up(i,j,k)+3.e0*up(i,j,k-1)-up(i,j,k-2))

            b=dzt22*(up(i,j,k+1)-2.e0*up(i,j,k)+up(i,j,k-1))

            c=dzt16*(3.e0*up(i,j,k)                                     &
     &        +2.e0*up(i,j,k+1)+up(i,j,k-2)-6.e0*up(i,j,k-1))

            uf(i,j,k)=up(i,j,k)+((a*wc8u+b)*wc8u+c)*wc8u

          else

            a=dzt36*(up(i,j,k+2)                                        &
     &        -3.e0*up(i,j,k+1)+3.e0*up(i,j,k)-up(i,j,k-1))

            b=dzt22*(up(i,j,k+1)-2.e0*up(i,j,k)+up(i,j,k-1))

            c=dzt16*(6.e0*up(i,j,k+1)                                   &
     &        -up(i,j,k+2)-2.e0*up(i,j,k-1)-3.e0*up(i,j,k))

            uf(i,j,k)=up(i,j,k)+((a*wc8u+b)*wc8u+c)*wc8u

          end if

        end do
        end do

!$omp end do

      end do

!$omp do schedule(runtime) private(i,j,wc8u,a,b,c)

      do j=2,nj-2
      do i=2,ni-1

        wc8u=.5e0*(wc8s(i-1,j,2)+wc8s(i,j,2))

        if(wc8u.gt.0.e0) then

          uf(i,j,2)=up(i,j,2)+dzt1*wc8u*(up(i,j,2)-up(i,j,1))

        else

          a=dzt36*(up(i,j,4)                                            &
     &      -3.e0*up(i,j,3)+3.e0*up(i,j,2)-up(i,j,1))

          b=dzt22*(up(i,j,3)-2.e0*up(i,j,2)+up(i,j,1))

          c=dzt16*(6.e0*up(i,j,3)                                       &
     &      -up(i,j,4)-2.e0*up(i,j,1)-3.e0*up(i,j,2))

          uf(i,j,2)=up(i,j,2)+((a*wc8u+b)*wc8u+c)*wc8u

        end if

        wc8u=.5e0*(wc8s(i-1,j,nkm2)+wc8s(i,j,nkm2))

        if(wc8u.gt.0.e0) then

          a=dzt36*(up(i,j,nkm1)                                         &
     &      -3.e0*up(i,j,nkm2)+3.e0*up(i,j,nkm3)-up(i,j,nkm4))

          b=dzt22*(up(i,j,nkm1)-2.e0*up(i,j,nkm2)+up(i,j,nkm3))

          c=dzt16*(3.e0*up(i,j,nkm2)                                    &
     &      +2.e0*up(i,j,nkm1)+up(i,j,nkm4)-6.e0*up(i,j,nkm3))

          uf(i,j,nkm2)=up(i,j,nkm2)+((a*wc8u+b)*wc8u+c)*wc8u

        else

          uf(i,j,nkm2)=up(i,j,nkm2)                                     &
     &      +dzt1*wc8u*(up(i,j,nkm1)-up(i,j,nkm2))

        end if

      end do
      end do

!$omp end do

! -----

! Calculate the v advection vertically.

      do k=3,nk-3

!$omp do schedule(runtime) private(i,j,wc8v,a,b,c)

        do j=2,nj-1
        do i=2,ni-2

          wc8v=.5e0*(wc8s(i,j-1,k)+wc8s(i,j,k))

          if(wc8v.gt.0.e0) then

            a=dzt36*(vp(i,j,k+1)                                        &
     &        -3.e0*vp(i,j,k)+3.e0*vp(i,j,k-1)-vp(i,j,k-2))

            b=dzt22*(vp(i,j,k+1)-2.e0*vp(i,j,k)+vp(i,j,k-1))

            c=dzt16*(3.e0*vp(i,j,k)                                     &
     &        +2.e0*vp(i,j,k+1)+vp(i,j,k-2)-6.e0*vp(i,j,k-1))

            vf(i,j,k)=vp(i,j,k)+((a*wc8v+b)*wc8v+c)*wc8v

          else

            a=dzt36*(vp(i,j,k+2)                                        &
     &        -3.e0*vp(i,j,k+1)+3.e0*vp(i,j,k)-vp(i,j,k-1))

            b=dzt22*(vp(i,j,k+1)-2.e0*vp(i,j,k)+vp(i,j,k-1))

            c=dzt16*(6.e0*vp(i,j,k+1)                                   &
     &        -vp(i,j,k+2)-2.e0*vp(i,j,k-1)-3.e0*vp(i,j,k))

            vf(i,j,k)=vp(i,j,k)+((a*wc8v+b)*wc8v+c)*wc8v

          end if

        end do
        end do

!$omp end do

      end do

!$omp do schedule(runtime) private(i,j,wc8v,a,b,c)

      do j=2,nj-1
      do i=2,ni-2

        wc8v=.5e0*(wc8s(i,j-1,2)+wc8s(i,j,2))

        if(wc8v.gt.0.e0) then

          vf(i,j,2)=vp(i,j,2)+dzt1*wc8v*(vp(i,j,2)-vp(i,j,1))

        else

          a=dzt36*(vp(i,j,4)                                            &
     &      -3.e0*vp(i,j,3)+3.e0*vp(i,j,2)-vp(i,j,1))

          b=dzt22*(vp(i,j,3)-2.e0*vp(i,j,2)+vp(i,j,1))

          c=dzt16*(6.e0*vp(i,j,3)                                       &
     &      -vp(i,j,4)-2.e0*vp(i,j,1)-3.e0*vp(i,j,2))

          vf(i,j,2)=vp(i,j,2)+((a*wc8v+b)*wc8v+c)*wc8v

        end if

        wc8v=.5e0*(wc8s(i,j-1,nkm2)+wc8s(i,j,nkm2))

        if(wc8v.gt.0.e0) then

          a=dzt36*(vp(i,j,nkm1)                                         &
     &      -3.e0*vp(i,j,nkm2)+3.e0*vp(i,j,nkm3)-vp(i,j,nkm4))

          b=dzt22*(vp(i,j,nkm1)-2.e0*vp(i,j,nkm2)+vp(i,j,nkm3))

          c=dzt16*(3.e0*vp(i,j,nkm2)                                    &
     &      +2.e0*vp(i,j,nkm1)+vp(i,j,nkm4)-6.e0*vp(i,j,nkm3))

          vf(i,j,nkm2)=vp(i,j,nkm2)+((a*wc8v+b)*wc8v+c)*wc8v

        else

          vf(i,j,nkm2)=vp(i,j,nkm2)                                     &
     &      +dzt1*wc8v*(vp(i,j,nkm1)-vp(i,j,nkm2))

        end if

      end do
      end do

!$omp end do

! -----

! Calculate the w advection vertically.

      do k=3,nk-2

!$omp do schedule(runtime) private(i,j,a,b,c)

        do j=2,nj-2
        do i=2,ni-2

          if(wc(i,j,k).gt.0.e0) then

            a=dzt36*(wp(i,j,k+1)                                        &
     &        -3.e0*wp(i,j,k)+3.e0*wp(i,j,k-1)-wp(i,j,k-2))

            b=dzt22*(wp(i,j,k+1)-2.e0*wp(i,j,k)+wp(i,j,k-1))

            c=dzt16*(3.e0*wp(i,j,k)                                     &
     &        +2.e0*wp(i,j,k+1)+wp(i,j,k-2)-6.e0*wp(i,j,k-1))

            wf(i,j,k)=wp(i,j,k)                                         &
     &        +((a*wc(i,j,k)+b)*wc(i,j,k)+c)*wc(i,j,k)

          else

            a=dzt36*(wp(i,j,k+2)                                        &
     &        -3.e0*wp(i,j,k+1)+3.e0*wp(i,j,k)-wp(i,j,k-1))

            b=dzt22*(wp(i,j,k+1)-2.e0*wp(i,j,k)+wp(i,j,k-1))

            c=dzt16*(6.e0*wp(i,j,k+1)                                   &
     &        -wp(i,j,k+2)-2.e0*wp(i,j,k-1)-3.e0*wp(i,j,k))

            wf(i,j,k)=wp(i,j,k)                                         &
     &        +((a*wc(i,j,k)+b)*wc(i,j,k)+c)*wc(i,j,k)

          end if

        end do
        end do

!$omp end do

      end do

!$omp do schedule(runtime) private(i,j,a,b,c)

      do j=2,nj-2
      do i=2,ni-2

        if(wc(i,j,2).gt.0.e0) then

          wf(i,j,2)=wp(i,j,2)+dzt1*wc(i,j,2)*(wp(i,j,2)-wp(i,j,1))

        else

          a=dzt36*(wp(i,j,4)                                            &
     &      -3.e0*wp(i,j,3)+3.e0*wp(i,j,2)-wp(i,j,1))

          b=dzt22*(wp(i,j,3)-2.e0*wp(i,j,2)+wp(i,j,1))

          c=dzt16*(6.e0*wp(i,j,3)                                       &
     &      -wp(i,j,4)-2.e0*wp(i,j,1)-3.e0*wp(i,j,2))

          wf(i,j,2)=wp(i,j,2)                                           &
     &      +((a*wc(i,j,2)+b)*wc(i,j,2)+c)*wc(i,j,2)

        end if

        if(wc(i,j,nkm1).gt.0.e0) then

          a=dzt36*(wp(i,j,nk)                                           &
     &      -3.e0*wp(i,j,nkm1)+3.e0*wp(i,j,nkm2)-wp(i,j,nkm3))

          b=dzt22*(wp(i,j,nk)-2.e0*wp(i,j,nkm1)+wp(i,j,nkm2))

          c=dzt16*(3.e0*wp(i,j,nkm1)                                    &
     &      +2.e0*wp(i,j,nk)+wp(i,j,nkm3)-6.e0*wp(i,j,nkm2))

          wf(i,j,nkm1)=wp(i,j,nkm1)                                     &
     &      +((a*wc(i,j,nkm1)+b)*wc(i,j,nkm1)+c)*wc(i,j,nkm1)

        else

          wf(i,j,nkm1)=wp(i,j,nkm1)                                     &
     &      +dzt1*wc(i,j,nkm1)*(wp(i,j,nk)-wp(i,j,nkm1))

        end if

      end do
      end do

!$omp end do

! -----

!$omp end parallel

!! -----

      end subroutine s_vculuvw

!-----7--------------------------------------------------------------7--

      end module m_vculuvw
