!***********************************************************************
      module m_vculs0
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/04/03
!     Modification: 2006/05/12, 2006/06/21, 2007/01/05, 2007/07/30,
!                   2007/10/19, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2009/02/27, 2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the vertical scalar advection by the Cubic Lagrange
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

      public :: vculs0, s_vculs0

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vculs0

        module procedure s_vculs0

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_vculs0(fpdziv,ivstp,dtsep,ni,nj,nk,wc8s,sp,sf)
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

      real, intent(in) :: wc8s(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity
                       ! at scalar points

! Input and output variables

      real, intent(inout) :: sp(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at past

      real, intent(inout) :: sf(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at future

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

        call copy3d(0,ni+1,0,nj+1,1,nk,sf,sp)

      end if

! -----

! Calculate the scalar advection vertically.

!$omp parallel default(shared) private(k)

      do k=3,nk-3

!$omp do schedule(runtime) private(i,j,a,b,c)

        do j=2,nj-2
        do i=2,ni-2

          if(wc8s(i,j,k).gt.0.e0) then

            a=dzt36*(sp(i,j,k+1)                                        &
     &        -3.e0*sp(i,j,k)+3.e0*sp(i,j,k-1)-sp(i,j,k-2))

            b=dzt22*(sp(i,j,k+1)-2.e0*sp(i,j,k)+sp(i,j,k-1))

            c=dzt16*(3.e0*sp(i,j,k)                                     &
     &        +2.e0*sp(i,j,k+1)+sp(i,j,k-2)-6.e0*sp(i,j,k-1))

            sf(i,j,k)=max(sp(i,j,k)                                     &
     &        +((a*wc8s(i,j,k)+b)*wc8s(i,j,k)+c)*wc8s(i,j,k),0.e0)

          else

            a=dzt36*(sp(i,j,k+2)                                        &
     &        -3.e0*sp(i,j,k+1)+3.e0*sp(i,j,k)-sp(i,j,k-1))

            b=dzt22*(sp(i,j,k+1)-2.e0*sp(i,j,k)+sp(i,j,k-1))

            c=dzt16*(6.e0*sp(i,j,k+1)                                   &
     &        -sp(i,j,k+2)-2.e0*sp(i,j,k-1)-3.e0*sp(i,j,k))

            sf(i,j,k)=max(sp(i,j,k)                                     &
     &        +((a*wc8s(i,j,k)+b)*wc8s(i,j,k)+c)*wc8s(i,j,k),0.e0)

          end if

        end do
        end do

!$omp end do

      end do

!$omp do schedule(runtime) private(i,j,a,b,c)

      do j=2,nj-2
      do i=2,ni-2

        if(wc8s(i,j,2).gt.0.e0) then

          sf(i,j,2)=max(sp(i,j,2)                                       &
     &      +dzt1*wc8s(i,j,2)*(sp(i,j,2)-sp(i,j,1)),0.e0)

        else

          a=dzt36*(sp(i,j,4)                                            &
     &      -3.e0*sp(i,j,3)+3.e0*sp(i,j,2)-sp(i,j,1))

          b=dzt22*(sp(i,j,3)-2.e0*sp(i,j,2)+sp(i,j,1))

          c=dzt16*(6.e0*sp(i,j,3)                                       &
     &      -sp(i,j,4)-2.e0*sp(i,j,1)-3.e0*sp(i,j,2))

          sf(i,j,2)=max(sp(i,j,2)                                       &
     &      +((a*wc8s(i,j,2)+b)*wc8s(i,j,2)+c)*wc8s(i,j,2),0.e0)

        end if

        if(wc8s(i,j,nkm2).gt.0.e0) then

          a=dzt36*(sp(i,j,nkm1)                                         &
     &      -3.e0*sp(i,j,nkm2)+3.e0*sp(i,j,nkm3)-sp(i,j,nkm4))

          b=dzt22*(sp(i,j,nkm1)-2.e0*sp(i,j,nkm2)+sp(i,j,nkm3))

          c=dzt16*(3.e0*sp(i,j,nkm2)                                    &
     &      +2.e0*sp(i,j,nkm1)+sp(i,j,nkm4)-6.e0*sp(i,j,nkm3))

          sf(i,j,nkm2)=max(sp(i,j,nkm2)                                 &
     &     +((a*wc8s(i,j,nkm2)+b)*wc8s(i,j,nkm2)+c)*wc8s(i,j,nkm2),0.e0)

        else

          sf(i,j,nkm2)=max(sp(i,j,nkm2)                                 &
     &      +dzt1*wc8s(i,j,nkm2)*(sp(i,j,nkm1)-sp(i,j,nkm2)),0.e0)

        end if

      end do
      end do

!$omp end do

!$omp end parallel

! -----

      end subroutine s_vculs0

!-----7--------------------------------------------------------------7--

      end module m_vculs0
