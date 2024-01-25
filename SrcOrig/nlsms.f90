!***********************************************************************
      module m_nlsms
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/07/15
!     Modification: 2003/11/28, 2003/12/12, 2004/02/01, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the non linear smoothing for optional scalar variable.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: nlsms, s_nlsms

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface nlsms

        module procedure s_nlsms

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
      subroutine s_nlsms(fpnlhcoe,fpnlvcoe,amp,ni,nj,nk,rbr,s,sfrc,rbrs,&
     &                   tmp1,tmp2,tmp3)
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
                       ! Amplitude of optional scalar variable

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: s(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable

! Input and output variable

      real, intent(inout) :: sfrc(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar forcing term

! Internal shared variables

      real nlhcoe      ! Horizontal non linear smoothig coefficient
      real nlvcoe      ! Vertical non linear smoothig coefficient

      real ampnlh      ! nlhcoe / amp
      real ampnlv      ! nlvcoe / amp

      real, intent(inout) :: rbrs(0:ni+1,0:nj+1,1:nk)
                       ! rbr x optional slcar perturbation

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

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpnlhcoe,nlhcoe)
      call getrname(fpnlvcoe,nlvcoe)

! -----

! Set the common used variables.

      ampnlh=nlhcoe/amp
      ampnlv=nlvcoe/amp

! -----

! Calculate the non linear scalar numerical smoothing.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          rbrs(i,j,k)=rbr(i,j,k)*s(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-2
        do i=2,ni-1
          a=rbrs(i,j,k)-rbrs(i-1,j,k)

          tmp1(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-1
        do i=2,ni-2
          a=rbrs(i,j,k)-rbrs(i,j-1,k)

          tmp2(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-2
        do i=2,ni-2
          a=rbrs(i,j,k)-rbrs(i,j,k-1)

          tmp3(i,j,k)=a*abs(a)

        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          sfrc(i,j,k)=sfrc(i,j,k)+(ampnlv*(tmp3(i,j,k+1)-tmp3(i,j,k))   &
     &      +ampnlh*((tmp1(i+1,j,k)-tmp1(i,j,k))                        &
     &      +(tmp2(i,j+1,k)-tmp2(i,j,k))))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_nlsms

!-----7--------------------------------------------------------------7--

      end module m_nlsms
