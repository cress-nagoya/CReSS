!***********************************************************************
      module m_lsps0
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/05/31
!     Modification: 2004/08/01, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the lateral sponge damping for optional scalar variable
!     to initial.

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

      public :: lsps0, s_lsps0

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface lsps0

        module procedure s_lsps0

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
      subroutine s_lsps0(fplspopt,fpwdnews,fplsnews,fplspsmt,ni,nj,nk,  &
     &                   rst,sp,rbcxy,sfrc,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpwdnews
                       ! Formal parameter of unique index of wdnews

      integer, intent(in) :: fplsnews
                       ! Formal parameter of unique index of lsnews

      integer, intent(in) :: fplspsmt
                       ! Formal parameter of unique index of lspsmt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jabobian

      real, intent(in) :: sp(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at past

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

! Input and output variable

      real, intent(inout) :: sfrc(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar forcing term

! Internal shared variables

      integer lspopt   ! Option for lateral sponge damping

      integer wdnews   ! Lateral sponge damping thickness

      real lsnews      ! Lateral sponge damping coefficient
      real lspsmt      ! Lateral sponge smoothing coefficient

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fplspopt,lspopt)
      call getiname(fpwdnews,wdnews)
      call getrname(fplsnews,lsnews)
      call getrname(fplspsmt,lspsmt)

! -----

! Calculate the lateral sponge damping for optional scalar variable to
! initial.

      if(wdnews.ge.1) then

!$omp parallel default(shared) private(k)

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            tmp1(i,j,k)=rst(i,j,k)*sp(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(lspopt.lt.10) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              sfrc(i,j,k)=sfrc(i,j,k)-lsnews*rbcxy(i,j)*tmp1(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j,a)

            do j=2,nj-2
            do i=2,ni-2
              a=2.e0*tmp1(i,j,k)

              sfrc(i,j,k)=sfrc(i,j,k)-rbcxy(i,j)*(lsnews*tmp1(i,j,k)    &
     &          -lspsmt*(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)              &
     &          +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a)))

            end do
            end do

!$omp end do

          end do

        end if

!$omp end parallel

      end if

! -----

      end subroutine s_lsps0

!-----7--------------------------------------------------------------7--

      end module m_lsps0
