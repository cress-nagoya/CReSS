!***********************************************************************
      module m_smoo2uvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/09/16
!     Modification: 1999/10/12, 1999/11/01, 1999/11/19, 1999/11/24,
!                   2000/01/17, 2001/06/06, 2002/04/02, 2002/08/15,
!                   2002/12/11, 2003/01/04, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/11/28, 2003/12/12, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the 2nd order smoothing for the velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: smoo2uvw, s_smoo2uvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface smoo2uvw

        module procedure s_smoo2uvw

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
      subroutine s_smoo2uvw(fpsmhcoe,fpsmvcoe,ni,nj,nk,jcb8w,           &
     &                      ubr,vbr,rbr,rst8w,u,v,w,ufrc,vfrc,wfrc,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsmhcoe
                       ! Formal parameter of unique index of smhcoe

      integer, intent(in) :: fpsmvcoe
                       ! Formal parameter of unique index of smvcoe

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

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

      real smhcoe      ! Horizontal smoothig coefficient
      real smvcoe      ! Vertical smoothig coefficient

      real smhc5       ! 0.5 x smhcoe
      real smvc5       ! 0.5 x smvcoe

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpsmhcoe,smhcoe)
      call getrname(fpsmvcoe,smvcoe)

! -----

! Set the common used variables.

      smhc5=.5e0*smhcoe
      smvc5=.5e0*smvcoe

! -----

!! Calculate the 2nd order velocity numerical smoothing.

!$omp parallel default(shared) private(k)

! Calculate the 2nd order u smoothing.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni
          tmp1(i,j,k)=(rbr(i-1,j,k)+rbr(i,j,k))*(u(i,j,k)-ubr(i,j,k))
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-2
        do i=2,ni-1
          a=2.e0*tmp1(i,j,k)

          ufrc(i,j,k)=ufrc(i,j,k)                                       &
     &      +(smvc5*((tmp1(i,j,k+1)+tmp1(i,j,k-1))-a)                   &
     &      +smhc5*(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)                   &
     &      +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a)))

        end do
        end do

!$omp end do

      end do

! -----

! Calculate the 2nd order v smoothing.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj
        do i=1,ni-1
          tmp1(i,j,k)=(rbr(i,j-1,k)+rbr(i,j,k))*(v(i,j,k)-vbr(i,j,k))
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-1
        do i=2,ni-2
          a=2.e0*tmp1(i,j,k)

          vfrc(i,j,k)=vfrc(i,j,k)                                       &
     &      +(smvc5*((tmp1(i,j,k+1)+tmp1(i,j,k-1))-a)                   &
     &      +smhc5*(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)                   &
     &      +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a)))

        end do
        end do

!$omp end do

      end do

! -----

! Calculate the 2nd order w smoothing.

      do k=1,nk

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          tmp1(i,j,k)=rst8w(i,j,k)*w(i,j,k)/jcb8w(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j,a)

        do j=2,nj-2
        do i=2,ni-2
          a=2.e0*tmp1(i,j,k)

          wfrc(i,j,k)=wfrc(i,j,k)                                       &
     &      +(smvcoe*((tmp1(i,j,k+1)+tmp1(i,j,k-1))-a)                  &
     &      +smhcoe*(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)                  &
     &      +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a)))

        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_smoo2uvw

!-----7--------------------------------------------------------------7--

      end module m_smoo2uvw
