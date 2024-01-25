!***********************************************************************
      module m_lspqv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/07
!     Modification: 1999/07/05, 1999/08/03, 1999/09/30, 1999/10/12,
!                   1999/11/01, 2000/01/17, 2000/02/07, 2000/04/18,
!                   2001/01/15, 2001/03/13, 2001/04/15, 2001/05/29,
!                   2001/06/06, 2001/06/29, 2001/07/13, 2001/08/07,
!                   2002/04/02, 2002/07/23, 2002/08/15, 2002/12/11,
!                   2003/01/04, 2003/04/30, 2003/05/19, 2003/06/27,
!                   2003/10/10, 2003/12/12, 2004/05/07, 2004/05/31,
!                   2004/08/01, 2004/08/20, 2006/09/21, 2007/05/07,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2009/03/23, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the lateral sponge damping for the water vapor mixing
!     ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: lspqv, s_lspqv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface lspqv

        module procedure s_lspqv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_lspqv(fpgpvvar,fplspopt,fpwdnews,fplsnews,fplspsmt,  &
     &                   gtinc,ni,nj,nk,rst,qvbr,qvp,rbcxy,qvgpv,qvtd,  &
     &                   qvfrc,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

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

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jabobian

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(in) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(in) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

! Input and output variable

      real, intent(inout) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

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

! Initialize the character variable.

      call inichar(gpvvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getiname(fplspopt,lspopt)
      call getiname(fpwdnews,wdnews)
      call getrname(fplsnews,lsnews)
      call getrname(fplspsmt,lspsmt)

! -----

!! Calculate the lateral sponge damping.

!$omp parallel default(shared) private(k)

      if(wdnews.ge.1) then

! Damp to the GPV data.

        if(gpvvar(2:2).eq.'o'.and.mod(lspopt,10).eq.1) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              tmp1(i,j,k)=rst(i,j,k)                                    &
     &          *(qvp(i,j,k)-(qvgpv(i,j,k)+qvtd(i,j,k)*gtinc))
            end do
            end do

!$omp end do

          end do

! -----

! Damp to the base state value.

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              tmp1(i,j,k)=rst(i,j,k)*(qvp(i,j,k)-qvbr(i,j,k))
            end do
            end do

!$omp end do

          end do

        end if

! -----

! Finally get the lateral sponge damping term.

        if(lspopt.lt.10) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              qvfrc(i,j,k)=qvfrc(i,j,k)-lsnews*rbcxy(i,j)*tmp1(i,j,k)
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

              qvfrc(i,j,k)=qvfrc(i,j,k)-rbcxy(i,j)*(lsnews*tmp1(i,j,k)  &
     &          -lspsmt*(((tmp1(i+1,j,k)+tmp1(i-1,j,k))-a)              &
     &          +((tmp1(i,j+1,k)+tmp1(i,j-1,k))-a)))

            end do
            end do

!$omp end do

          end do

        end if

! -----

      end if

!$omp end parallel

!! -----

      end subroutine s_lspqv

!-----7--------------------------------------------------------------7--

      end module m_lspqv
