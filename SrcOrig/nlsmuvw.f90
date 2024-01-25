!***********************************************************************
      module m_nlsmuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/07/15
!     Modification: 2003/11/28, 2003/12/12, 2004/02/01, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the non linear smoothing for the velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: nlsmuvw, s_nlsmuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface nlsmuvw

        module procedure s_nlsmuvw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_nlsmuvw(fpnlhcoe,fpnlvcoe,amp,ni,nj,nk,jcb8w,        &
     &                     ubr,vbr,rbr,rst8w,u,v,w,ufrc,vfrc,wfrc,      &
     &                     tmp1,tmp2,tmp3,tmp4)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpnlhcoe
                       ! Formal parameter of unique index of nlhcoe

      integer, intent(in) :: fpnlvcoe
                       ! Formal parameter of unique index of nlvcoe

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: amp
                       ! Amplitude of velocity

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

! Input and output variables

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      real nlhcoe      ! Horizontal non linear smoothig coefficient
      real nlvcoe      ! Vertical non linear smoothig coefficient

      real ampnlh      ! nlhcoe / amp
      real ampnlv      ! nlvcoe / amp

      real amph25      ! 0.25 x nlhcoe / amp
      real ampv25      ! 0.25 x nlvcoe / amp

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpnlhcoe,nlhcoe)
      call getrname(fpnlvcoe,nlvcoe)

! -----

! Set the common used variables.

      ampnlh=nlhcoe/amp
      ampnlv=nlvcoe/amp

      amph25=.25e0*nlhcoe/amp
      ampv25=.25e0*nlvcoe/amp

! -----

!! Calculate the non linear velocity numerical smoothing.

!$omp parallel default(shared) private(k)

! Calculate the non linear u smoothing.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni
          tmp4(i,j,k)=(rbr(i-1,j,k)+rbr(i,j,k))*(u(i,j,k)-ubr(i,j,k))
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-2
        do i=1,ni-1
          a=tmp4(i+1,j,k)-tmp4(i,j,k)

          tmp1(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-1
        do i=2,ni-1
          a=tmp4(i,j,k)-tmp4(i,j-1,k)

          tmp2(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-2
        do i=2,ni-1
          a=tmp4(i,j,k)-tmp4(i,j,k-1)

          tmp3(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          ufrc(i,j,k)=ufrc(i,j,k)+(ampv25*(tmp3(i,j,k+1)-tmp3(i,j,k))   &
     &      +amph25*((tmp1(i,j,k)-tmp1(i-1,j,k))                        &
     &      +(tmp2(i,j+1,k)-tmp2(i,j,k))))
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the non linear v smoothing.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj
        do i=1,ni-1
          tmp4(i,j,k)=(rbr(i,j-1,k)+rbr(i,j,k))*(v(i,j,k)-vbr(i,j,k))
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-1
        do i=2,ni-1
          a=tmp4(i,j,k)-tmp4(i-1,j,k)

          tmp1(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j,a)

        do j=1,nj-1
        do i=2,ni-2
          a=tmp4(i,j+1,k)-tmp4(i,j,k)

          tmp2(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-1
        do i=2,ni-2
          a=tmp4(i,j,k)-tmp4(i,j,k-1)

          tmp3(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          vfrc(i,j,k)=vfrc(i,j,k)+(ampv25*(tmp3(i,j,k+1)-tmp3(i,j,k))   &
     &      +amph25*((tmp1(i+1,j,k)-tmp1(i,j,k))                        &
     &      +(tmp2(i,j,k)-tmp2(i,j-1,k))))
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the non linear w smoothing.

      do k=1,nk

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          tmp4(i,j,k)=rst8w(i,j,k)*w(i,j,k)/jcb8w(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-2
        do i=2,ni-1
          a=tmp4(i,j,k)-tmp4(i-1,j,k)

          tmp1(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-1
        do i=2,ni-2
          a=tmp4(i,j,k)-tmp4(i,j-1,k)

          tmp2(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

      end do

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-2
        do i=2,ni-2
          a=tmp4(i,j,k+1)-tmp4(i,j,k)

          tmp3(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          wfrc(i,j,k)=wfrc(i,j,k)+(ampnlv*(tmp3(i,j,k)-tmp3(i,j,k-1))   &
     &      +ampnlh*((tmp1(i+1,j,k)-tmp1(i,j,k))                        &
     &      +(tmp2(i,j+1,k)-tmp2(i,j,k))))
        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_nlsmuvw

!-----7--------------------------------------------------------------7--

      end module m_nlsmuvw
